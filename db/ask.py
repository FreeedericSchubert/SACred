"""Ask the database in plain language: Claude writes DuckDB SQL, runs it read-only on the release tables, answers.

    python db/ask.py "Which graphene sites have a DFT dG_H within 0.2 eV of zero?"
    python db/ask.py --prompt          # print the system prompt (instructions and the schema Claude gets)

The explorer (notebooks/sac_explorer.py) puts the same agent behind its chat box. Claude (MODEL) gets the schema
from column_docs and one tool, run_sql. The tables are views on the release Parquet files (data/release/parquet)
in an in-memory DuckDB that may read exactly those files and nothing else, and only a single SELECT gets through.
Needs an Anthropic API key: ANTHROPIC_API_KEY in the environment, or for the explorer a line
ANTHROPIC_API_KEY=... in notebooks/.env (marimo loads that file at start; .env is git-ignored).
"""
import asyncio
import numbers
import os
import sys
import time

import anthropic
import duckdb
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REL = os.path.join(ROOT, "data", "release", "parquet")
MODEL = "claude-opus-5-5"
EFFORT = "low"                 # short turns for a live demo; "medium" thinks longer on hard questions
FALLBACK_BETA = "server-side-fallback-2026-07-01"   # a declined request is re-run on the recommended fallback model
MAX_TOOL_TURNS = 8             # rounds of queries per question
QUERY_SECONDS = 20             # a query running longer is interrupted
ROWS_KEPT = 200                # rows fetched per query
ROWS_TO_MODEL = 50             # rows Claude sees (CSV)
ROWS_SHOWN = 8                 # rows shown in the chat

EXAMPLES = [
    "Which graphene SACs have a DFT ΔG_H within ±0.2 eV, and which enzyme sites are they modelled on?",
    "How many unique sites does each enzyme metal have, and how many come from HER, ORR, CO₂RR or NRR enzymes?",
    "Compare Co and Ni on graphene, MoS₂ and silica by xTB ΔG_H.",
    "Which enzyme first shells fit MoS₂ best?",
]

NO_KEY = ("**No Anthropic API key found.** Put a line `ANTHROPIC_API_KEY=...` (a key from the Claude Console) into "
          "`notebooks/.env`, which git ignores, and restart the explorer. The SQL box below works without it.")

SYSTEM = """You answer questions about a research database of metal sites in proteins and of the single-atom \
catalysts (SACs) derived from them. You get the numbers by querying the database with DuckDB SQL through the run_sql \
tool.

The project: every d-block metal site in the PDB was extracted (sites), its first coordination shell described \
(donor atoms, geometry) and mapped onto a single-atom site on a support: a double (DV, four in-plane donors) or single \
(SV, three donors) vacancy in graphene with the donor atoms as dopants (support_sites has the other supports: a metal \
on three S of MoS2, a metal in an S vacancy of MoS2 (MoS2vS), a metal on siloxide O of silica (SiO2)). For a \
selection of these SAC models the hydrogen adsorption free energy dG_H was computed with GFN2-xTB (a fast screen: \
xtb_scan for graphene, xtb_supports for MoS2 and silica) and with DFT (Quantum ESPRESSO, PBE-D3; dft_results). A \
literature check (literature_check, literature_papers) tested whether a motif had been reported before.

Conventions:
- dG_H is in eV, dG_H = dE_H + 0.24 eV. Near 0 is ideal for hydrogen evolution (HER); |dG_H| <= 0.2 eV (DFT) counts \
as promising, and the explorer uses |dG_H| <= 0.3 eV for the coarser xTB values. Negative: H binds too strongly; \
positive: too weakly. DFT is the reference. xTB values come raw (dG_H) and with a per-metal offset fitted on the DFT \
pairs (xtb_scan.dG_H_calibrated, xtb_supports.dG_H_offset). Always say which kind of value you report.
- state: H_M = H on the metal, H_S = H on an S donor, H_O = H on a siloxide O. A model can have several states; \
report them separately or say which one you use.
- sites has one row per metal atom. Copies of the same site in one entry share a copy_group; count unique sites \
with WHERE is_representative. metal_class = 'biological' are the nine enzyme metals V, Mn, Fe, Co, Ni, Cu, Zn, Mo, \
W; the other metals are mostly soaks, drugs and phasing agents.
- site_sets.set_name: 'complete' (all sites), 'enzyme' (sites of HER, ORR, CO2RR and NRR enzymes), 'pdb30' (one \
entry per 30 % sequence cluster). sites.reaction holds those labels ';'-joined (e.g. 'CO2RR;HER') and '' when \
there is none, so filter with reaction LIKE '%HER%'.
- Motif notation: sites.graphene_motif reads like 'Ni-S2C2(DV,cis)' or 'Os-N4(DV)+ax:N2' (axial donors after \
+ax:). The SAC models (sac_models, xtb_scan, dft_results, dft_stability) use model_id 'graphene:Ni-N2S2:DV:cis5' and \
motif 'Ni-N2S2:DV:cis5'. Link a model to the enzyme sites it stands for with sac_models.sweep_motif = \
split_part(sites.graphene_motif, '+', 1); sac_models.origin names the enzyme site a model was built from.
- support_motifs summarises the mapping onto MoS2, MoS2vS and SiO2 per motif (frac_fit: share of sites whose donor \
geometry fits the support).

How to work:
- Every number in your answer comes from run_sql; never invent values. If you are unsure about a table, DESCRIBE it \
or look at a few rows first.
- Keep queries small: aggregate, filter, and LIMIT 20 unless the question needs more rows. Only a single SELECT runs \
(WITH, DESCRIBE and SUMMARIZE included).
- If a query fails, read the error, fix the query and run it again.
- If the question is ambiguous, take the most sensible reading and say in one line which one you took.
- The chat shows every query you run and the first rows of its result, so do not repeat the SQL in the answer.
- Answer in the language of the question: a few sentences of Markdown, with a small table (at most 10 rows) only \
when it helps. Round energies to 0.01 eV and distances to 0.01 Å, and give units.
- If the database cannot answer the question, say so and say what it does contain.

Tables (DuckDB views: name, rows, then each column with its type and meaning):
{schema}
"""

TOOLS = [{
    "name": "run_sql",
    "description": "Run one read-only DuckDB SQL query on the SAC database and get the result as CSV (up to "
                   f"{ROWS_TO_MODEL} rows, with the number of rows). Use it for every number, count, ranking or "
                   "example in your answer: the answer must come from the data. Only a single SELECT runs (WITH, "
                   "DESCRIBE and SUMMARIZE included). If it returns an error, fix the query and call it again.",
    "input_schema": {"type": "object", "additionalProperties": False, "required": ["sql"],
                     "properties": {"sql": {"type": "string", "description": "one DuckDB SELECT query"}}},
    "strict": True,
}]


def has_credentials():
    """Whether the SDK will find credentials (checked without reading them): a key or token in the environment, or
    a profile from `ant auth login`."""
    return bool(os.environ.get("ANTHROPIC_API_KEY") or os.environ.get("ANTHROPIC_AUTH_TOKEN")
                or os.path.isdir(os.path.expanduser("~/.config/anthropic")))


# ---------------------------------------------------------------- database

class Database:
    """The release tables as views. DuckDB may read exactly these Parquet files, the settings are locked, and only a
    single SELECT is passed on (COPY would otherwise be able to write into the release folder)."""

    def __init__(self, rel=REL):
        self.con = duckdb.connect()
        files = sorted(f for f in os.listdir(rel) if f.endswith(".parquet"))
        self.tables = [f[:-8] for f in files]
        paths = [os.path.join(rel, f) for f in files]
        for table, path in zip(self.tables, paths):
            self.con.execute(f"CREATE VIEW {table} AS SELECT * FROM read_parquet('{path}')")
        self.con.execute("SET allowed_paths = [%s]" % ", ".join("'%s'" % p for p in paths))
        self.con.execute("SET threads = 4")
        self.con.execute("SET memory_limit = '2GB'")
        self.con.execute("SET enable_external_access = false")
        self.con.execute("SET lock_configuration = true")

    def schema(self):
        """Every table with its row count and every column with its type and its line from column_docs."""
        docs = self.con.execute("SELECT table_name, column_name, description FROM column_docs").fetchall()
        doc = {(t, c): d for t, c, d in docs}
        lines = []
        for table in self.tables:
            n = self.con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            lines.append(f"\n{table} ({n:,} rows)")
            for column, kind, *_ in self.con.execute(f"DESCRIBE {table}").fetchall():
                text = doc.get((table, column))
                lines.append(f"  {column} {kind}" + (f": {text}" if text else ""))
        return "\n".join(lines).strip()

    async def query(self, sql):
        """-> (DataFrame of at most ROWS_KEPT rows, whether there were more). Raises on anything but one SELECT, on
        DuckDB errors and after QUERY_SECONDS."""
        statements = duckdb.extract_statements(sql)
        if len(statements) != 1 or statements[0].type != duckdb.StatementType.SELECT:
            raise ValueError("only a single SELECT query runs here (WITH, DESCRIBE and SUMMARIZE included)")
        cursor = self.con.cursor()
        try:
            return await asyncio.wait_for(asyncio.to_thread(self._fetch, cursor, sql), QUERY_SECONDS)
        except asyncio.TimeoutError:
            cursor.interrupt()
            raise TimeoutError(f"the query ran longer than {QUERY_SECONDS} s and was stopped; narrow it down")

    @staticmethod
    def _fetch(cursor, sql):
        cursor.execute(sql)
        columns = [d[0] for d in cursor.description]
        rows = cursor.fetchmany(ROWS_KEPT + 1)
        return pd.DataFrame(rows[:ROWS_KEPT], columns=columns), len(rows) > ROWS_KEPT


# ---------------------------------------------------------------- formatting

def _rows(df, more):
    return f"{len(df)}{'+' if more else ''} row{'' if len(df) == 1 and not more else 's'}"


def for_model(df, more):
    """The result as Claude sees it: row count, then CSV of the first ROWS_TO_MODEL rows."""
    head = _rows(df, more) + (f", the first {ROWS_TO_MODEL} below" if len(df) > ROWS_TO_MODEL else "")
    return (head + "\n" + df.head(ROWS_TO_MODEL).to_csv(index=False, float_format="%.4g"))[:12000]


def _cell(v):
    if v is None or (isinstance(v, float) and v != v):
        return ""
    if isinstance(v, (bool, str)) or not isinstance(v, numbers.Number):
        s = str(v).replace("|", "\\|").replace("\n", " ")
        return s if len(s) <= 48 else s[:47] + "…"
    if isinstance(v, numbers.Integral):
        return f"{v:,}"
    return f"{v:,.0f}" if abs(v) >= 1000 else f"{v:.3g}"


def for_chat(df, more):
    """A Markdown preview of the result: the first ROWS_SHOWN rows and 8 columns."""
    if not len(df.columns) or df.empty:
        return "_no rows_"
    d = df.iloc[:ROWS_SHOWN, :8]
    lines = ["| " + " | ".join(_cell(c) for c in d.columns) + " |", "|" + "---|" * len(d.columns)]
    lines += ["| " + " | ".join(_cell(v) for v in row) + " |" for row in d.itertuples(index=False)]
    if len(df) > ROWS_SHOWN or more:
        lines.append(f"\n_{_rows(df, more)} in all_")
    return "\n".join(lines)


def echo(content):
    """The response content to send back. After a fallback, the declined attempt's thinking and tool calls (blocks
    before the last fallback block) are left out, its text kept; without a fallback everything goes back unchanged."""
    last = max((i for i, b in enumerate(content) if b.type == "fallback"), default=-1)
    return [b for i, b in enumerate(content) if i > last or b.type == "text"]


# ---------------------------------------------------------------- agent

class Agent:
    """One conversation about the database. The message history only grows (earlier turns are never edited), so
    Claude's thinking blocks stay valid and the prompt cache keeps hitting."""

    def __init__(self, database=None):
        self.db = database or Database()
        self.system = SYSTEM.replace("{schema}", self.db.schema())
        self.client = None
        self.fallbacks = True
        self.history = []          # messages as sent to the API
        self.asked = []            # the questions in it

    async def _create(self):
        extra = dict(betas=[FALLBACK_BETA], fallbacks="default") if self.fallbacks else {}
        return await self.client.beta.messages.create(
            model=MODEL, max_tokens=16000, thinking={"type": "adaptive"}, output_config={"effort": EFFORT},
            system=self.system, tools=TOOLS, messages=self.history, cache_control={"type": "ephemeral"}, **extra)

    async def ask(self, question):
        """Async generator of Markdown pieces: each query with a preview of its result, then the answer."""
        if not has_credentials():
            yield NO_KEY
            return
        if self.client is None:
            self.client = anthropic.AsyncAnthropic()
        self.history.append({"role": "user", "content": question})
        self.asked.append(question)
        start, n = time.time(), 0
        for _ in range(MAX_TOOL_TURNS + 1):
            try:
                try:
                    r = await self._create()
                except anthropic.BadRequestError as e:
                    if not (self.fallbacks and "fallback" in str(e).lower()):
                        raise
                    self.fallbacks = False         # the account does not take the fallback beta: ask without it
                    r = await self._create()
            except TypeError as e:                 # the SDK found no credentials after all
                yield NO_KEY if "authentication" in str(e) else f"\n\n**Error:** {e}"
                return
            except anthropic.AuthenticationError:
                yield "\n\n**The API key was not accepted** (401). Check ANTHROPIC_API_KEY in notebooks/.env."
                return
            except anthropic.PermissionDeniedError as e:
                yield f"\n\n**The key has no access** to {MODEL} or a beta it needs: {e.message}"
                return
            except anthropic.RateLimitError:
                yield "\n\n**Rate limited** by the API; ask again in a minute."
                return
            except anthropic.APIStatusError as e:
                yield f"\n\n**The API answered {e.status_code}:** {e.message}"
                return
            except anthropic.APIConnectionError:
                yield "\n\n**No connection to the Anthropic API.** Check the network and ask again."
                return
            if r.stop_reason == "refusal":         # declined (also by the fallback): the turn is not kept
                yield "\n\nThe model declined this question; try asking it differently."
                return
            calls = [b for b in r.content if b.type == "tool_use"]
            notes = "\n\n".join(b.text.strip() for b in r.content if b.type == "text" and b.text.strip())
            if r.stop_reason == "max_tokens" and calls:   # cut off inside a tool call: not kept
                yield "\n\nThe answer was cut off; ask a narrower question."
                return
            self.history.append({"role": "assistant", "content": echo(r.content)})
            if r.stop_reason == "pause_turn":
                continue
            if r.stop_reason != "tool_use" or not calls:
                by = "" if r.model == MODEL else f" · answered by {r.model}"
                took = f"{n} quer{'y' if n == 1 else 'ies'} · {time.time() - start:.0f} s{by}"
                yield f"\n\n{notes}\n\n_{took}_"
                return
            if notes:
                yield f"\n\n{notes}"
            results = []
            for call in calls:
                sql = call.input.get("sql", "") if isinstance(call.input, dict) else ""
                n += 1
                try:
                    df, more = await self.db.query(sql)
                    content, shown, failed = for_model(df, more), for_chat(df, more), False
                    head = f"**Query {n}** · {_rows(df, more)}"
                except Exception as e:             # a wrong query: Claude reads the error and tries again
                    content, shown, failed = f"Error: {e}", "", True
                    head = f"**Query {n}** · failed: {str(e).splitlines()[0][:160]}"
                yield f"\n\n{head}\n\n```sql\n{sql.strip()}\n```\n\n{shown}"
                results.append({"type": "tool_result", "tool_use_id": call.id, "content": content, "is_error": failed})
            self.history.append({"role": "user", "content": results})
        yield f"\n\nStopped after {n} queries without an answer; ask a narrower question."

    async def chat(self, messages, config=None):
        """The model function of mo.ui.chat: answers the last message; marimo passes the whole chat each time."""
        asked = [str(m.content) for m in messages if m.role == "user"]
        if asked[:-1] != self.asked:               # the chat was cleared or edited: continue from what it shows
            self.history = [{"role": m.role, "content": str(m.content)} for m in messages[:-1]
                            if m.role in ("user", "assistant") and str(m.content).strip()]
            self.asked = asked[:-1]
        async for piece in self.ask(asked[-1]):
            yield piece


async def _main(question):
    async for piece in Agent().ask(question):
        print(piece, end="", flush=True)
    print()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    if sys.argv[1] == "--prompt":
        print(Agent().system)
    else:
        asyncio.run(_main(" ".join(sys.argv[1:])))
