"""Verify source/text hashes only; never builds or runs MHD."""
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
m=json.loads((root/'SOURCE_MANIFEST.json').read_text())
for path,r in {**m['source_files'],**m['header_files']}.items():
 assert hashlib.sha256((root/path).read_bytes()).hexdigest()==r['sha256'],path
print('PASS:',len(m['source_files']),'executed translation units;',len(m['header_files']),'headers; accepted full12h',m['whole_process_RT'],'xRT; portable rebuild remains unqualified')
