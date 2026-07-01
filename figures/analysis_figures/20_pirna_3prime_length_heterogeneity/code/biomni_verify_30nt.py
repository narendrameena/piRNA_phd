#!/usr/bin/env python3
"""DOUBLE-VERIFY across BioMNI agents: which protein produces/cuts the ~30 nt PACHYTENE piRNA class in MOUSE?"""
import asyncio
from fastmcp import Client
Q=("In MOUSE (Mus musculus) postnatal/adult testis, DOUBLE-VERIFY the biology of the ~30 nt PACHYTENE piRNA class "
   "specifically. CONFIRM or CORRECT each, with mechanism and key mouse references:\n"
   "(1) which PIWI protein BINDS the ~29-30 nt pachytene piRNAs? (expected: MIWI / PIWIL1)\n"
   "(2) which nuclease performs the 3'->5' TRIMMING that sets the mature ~30 nt 3' end to the bound-PIWI footprint? "
   "(expected: PNLDC1 / Trimmer, PARN-like)\n"
   "(3) what TERMINATES trimming / stabilizes the 3' end? (expected: HENMT1 2'-O-methylation)\n"
   "(4) what drives the POSTNATAL onset of the 30 nt pachytene class (~P14 onward)? (expected: A-MYB / MYBL1)\n"
   "(5) in Pnldc1-knockout mouse testis, are pachytene piRNAs LONGER/untrimmed (~31-35 nt) - confirming PNLDC1 "
   "sets the ~30 nt mature length, and that the 3' end is otherwise ragged/imprecise?\n"
   "For EACH: established IN MOUSE? mechanism? key references. Be specific to the 30 nt pachytene/MIWI class.")
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
        await asyncio.sleep(45)
asyncio.run(main())
