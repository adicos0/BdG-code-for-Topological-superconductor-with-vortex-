"""Verify provenance of the source-only notebook edition without running Wolfram."""
from pathlib import Path
import hashlib
import json
import re
import nbformat
from extract_original_code import extract

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((ROOT/'reference_pages/manifest.json').read_text())
    for name, expected in manifest['original_files_sha256'].items():
        assert sha256(ROOT/name) == expected, f'Original changed: {name}'
    sources = {}
    for stem in ['bdg_topological', 'added_vortex']:
        saved = json.loads((ROOT/'original_code'/f'{stem}.json').read_text())
        fresh = extract(ROOT/saved['source_file'])
        assert saved['cells'] == fresh, f'Transcription differs: {stem}'
        sources[stem] = fresh
        text = '\n\n'.join(f'(* Original input {c["index"]}; {saved["source_file"]}, line {c["line"]} *)\n'+c['code'] for c in fresh)+'\n'
        assert (ROOT/'original_code'/f'{stem}.wl').read_text() == text
    for page in manifest['pages']:
        assert sha256(ROOT/page['path']) == page['sha256'], page['path']
    code_count, page_count = 0, 0
    for path in sorted((ROOT/'notebooks').glob('*.ipynb')):
        notebook = nbformat.read(path, as_version=4)
        nbformat.validate(notebook)
        for cell in notebook.cells:
            assert cell.cell_type == 'markdown', 'Source-only edition has no execution outputs'
            if 'original_code' in cell.metadata:
                info = cell.metadata['original_code']
                original = sources[info['source_key']][info['input']-1]['code']
                expected = original[info['start']:info['end']]
                actual = cell.source.split('```mathematica\n', 1)[1].rsplit('\n```', 1)[0]
                assert actual == expected, f'Code changed: {path.name}'
                code_count += 1
            if 'thesis_source' in cell.metadata:
                info = cell.metadata['thesis_source']
                assert info['pdf_page'] == info['printed_page']+14
                assert (ROOT/info['asset']).is_file(), info['asset']
                page_count += 1
            for referenced in re.findall(r'/reference_pages/(thesis_p\d+\.png)', cell.source):
                assert (ROOT/'reference_pages'/referenced).is_file()
        print('PASS', path.name)
    print(f'{code_count} exact code excerpts; {page_count} page placements; {len(manifest["pages"])} unique page images.')
    print('Original .nb/.wdx hashes, BoxData transcription and notebook schema verified.')


if __name__ == '__main__':
    main()
