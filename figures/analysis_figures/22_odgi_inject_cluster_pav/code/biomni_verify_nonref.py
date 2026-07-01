#!/usr/bin/env python3
"""Triple-verify the non-reference / new piRNA-cluster biology across the 3 BioMNI agents."""
import asyncio
from fastmcp import Client
Q=("In mouse testis piRNA biology, please CONFIRM or CORRECT each statement, with the molecular mechanism and key "
   "references:\n"
   "(1) NEW / strain-specific piRNA clusters arise predominantly from transposable-element (TE) insertions (young "
   "LTR/ERVK, LINE-1, ERVL/MaLR) that occur AFTER mouse strains/lineages diverge.\n"
   "(2) a new active TE insertion can SEED a NEW piRNA cluster via the host's ADAPTIVE piRNA response (the piRNA-TE "
   "arms race / co-option), producing piRNAs that silence the element.\n"
   "(3) such young, strain-specific TE-insertion-driven clusters are ABSENT from the C57BL/6J (GRCm39) reference "
   "genome and can be WELL-EXPRESSED.\n"
   "(4) wild-derived / divergent strains (SPRET, CAST, PWK) carry MORE such non-reference clusters due to more recent "
   "TE activity.\nState for EACH whether established IN MOUSE, the mechanism, and key references.")
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
