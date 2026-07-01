#!/usr/bin/env python3
"""Triple-verify the piRNA length/3'-trimming biology across the 3 BioMNI agents (genomics/literature/general)."""
import asyncio
from fastmcp import Client
Q=("In mouse fetal and postnatal testis piRNA biogenesis, please CONFIRM or CORRECT each statement, with the "
   "molecular mechanism and key references:\n"
   "(1) piRNA length is set by the bound PIWI protein's footprint PLUS 3'->5' exonucleolytic trimming by PNLDC1 "
   "(Trimmer; PARN-like), with 2'-O-methylation by HENMT1 terminating trimming.\n"
   "(2) the loaded PIWI sets the mature length: MILI/PIWIL2 ~26-27 nt (pre-pachytene; fetal/perinatal), "
   "MIWI2/PIWIL4 ~28 nt (fetal, nuclear), MIWI/PIWIL1 ~29-30 nt (pachytene; postnatal/adult).\n"
   "(3) there is a developmental length SHIFT from ~26-27 nt (pre-pachytene, ~E16.5/perinatal) to ~29-30 nt "
   "(pachytene, ~P20.5/adult) across mouse spermatogenesis.\n"
   "(4) piRNA 3' ends are heterogeneous/ragged: a DEFINED 5' end but VARIABLE 3' length, because trimming is not "
   "single-nucleotide precise, so a shorter read is often the EXACT 5' prefix of a longer piRNA isoform.\n"
   "State for EACH whether it is established, with references.")
SRV=[("genomics",18881,"run_genomics_query"),("literature",18882,"run_literature_query"),("general",18883,"run_biomedical_query")]
async def one(name,port,tool):
    try:
        async with Client(f"http://127.0.0.1:{port}/mcp") as c:
            r=await c.call_tool(tool,{"query":Q,"use_cache":True})
            out=getattr(r,"data",None) or "\n".join(getattr(b,"text",str(b)) for b in (getattr(r,"content",None) or []))
            return name,out
    except Exception as e: return name,f"ERROR: {e}"
async def main():
    for s in SRV:
        name,out=await one(*s); print(f"\n{'='*22} BioMNI: {name} {'='*22}\n{out}",flush=True)
        await asyncio.sleep(45)   # sequential -> respect Groq 30k TPM (parallel tripped it)
asyncio.run(main())
