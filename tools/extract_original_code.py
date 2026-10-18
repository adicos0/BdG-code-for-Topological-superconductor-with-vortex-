"""Read original Wolfram RowBox input cells without evaluating or translating them.

Only strings and RowBox lists are accepted. Unknown box types fail explicitly.
Original BoxData is saved alongside the readable form for direct comparison.
"""
from pathlib import Path
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]


def closing_bracket(source, start):
    depth, quoted, escaped = 0, False, False
    for index in range(start, len(source)):
        char = source[index]
        if quoted:
            if escaped:
                escaped = False
            elif char == '\\':
                escaped = True
            elif char == '"':
                quoted = False
        elif char == '"':
            quoted = True
        elif char == '[':
            depth += 1
        elif char == ']':
            depth -= 1
            if depth == 0:
                return index
    raise ValueError('Unclosed box expression')


def readable_box(box):
    # Escaped newline continuations must remain inside a single quoted token.
    tokens = re.findall(r'"(?:\\[\s\S]|[^"\\])*"|[A-Za-z$][A-Za-z0-9$]*|[^\s]', box)
    position = 0

    def parse():
        nonlocal position
        token = tokens[position]
        position += 1
        if token.startswith('"'):
            value = token[1:-1].replace('\\\n', '')
            # Decode only the outer notebook string's escapes. Preserve named
            # Wolfram character escapes and inner plot-label formatting.
            return re.sub(r'\\(\\|"|n|r|t)',
                          lambda m: {'\\':'\\', '"':'"', 'n':'\n', 'r':'\r', 't':'\t'}[m[1]], value)
        if token == '{':
            values = []
            while tokens[position] != '}':
                values.append(parse())
                if tokens[position] == ',':
                    position += 1
            position += 1
            return values
        if token == 'RowBox':
            assert tokens[position] == '['
            position += 1
            values = parse()
            assert tokens[position] == ']'
            position += 1
            return ''.join(values)
        raise ValueError(f'Unsupported box token: {token!r}')

    result = parse()
    assert position == len(tokens), 'Trailing box syntax'
    # These named characters are whitespace in the displayed source input.
    return result.replace('\\[IndentingNewLine]', '\n')


def extract(path):
    source = path.read_text()
    cells = []
    for match in re.finditer(r'Cell\[BoxData\[', source):
        start = match.end()-1
        end = closing_bracket(source, start)
        if not re.match(r'\s*,\s*"(?:Input|Code)"', source[end+1:]):
            continue
        box = source[start+1:end]
        # The trailing empty GraphicsBox is a drawing canvas, not a program.
        if box.lstrip().startswith('GraphicsBox['):
            continue
        code = readable_box(box)
        if not code.strip():
            continue
        cells.append({'index': len(cells)+1,
                      'line': source[:match.start()].count('\n')+1,
                      'boxdata': box,
                      'code': code,
                      'code_sha256': hashlib.sha256(code.encode()).hexdigest()})
    return cells


def main():
    destination = ROOT/'original_code'
    destination.mkdir(exist_ok=True)
    for name in ['BdG topological.nb', 'added vortex.nb']:
        cells = extract(ROOT/name)
        stem = 'bdg_topological' if name.startswith('BdG') else 'added_vortex'
        (destination/f'{stem}.json').write_text(json.dumps({'source_file': name, 'cells': cells}, indent=2))
        text = '\n\n'.join(f'(* Original input {c["index"]}; {name}, line {c["line"]} *)\n'+c['code'] for c in cells)
        (destination/f'{stem}.wl').write_text(text+'\n')
        print(name, len(cells), 'nonempty calculation input cells')


if __name__ == '__main__':
    main()
