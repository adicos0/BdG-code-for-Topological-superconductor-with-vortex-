"""Extract only existing diagram rectangles, keeping every plot pixel intact.

The input PNGs are preserved page renders. No rescaling, filtering, recolouring,
annotations or numerical regeneration is applied to the diagram regions.
"""
from pathlib import Path
import json
import hashlib
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
# (asset name, printed source page, rectangle: left, top, right, bottom)
REGIONS = [
 ('lattice_connections',19,(195,110,825,600)),
 ('spectrum_zero_pairing',32,(190,115,835,562)),
 ('ldos_central_peak',33,(150,345,875,815)),
 ('spectrum_gap',34,(150,195,875,695)),
 ('ldos_gap',35,(150,110,875,575)),
 ('spectrum_in_gap',36,(155,390,875,889)),
 ('pairing_line_profile',37,(155,190,875,650)),
 ('pairing_map_smooth',38,(230,110,800,675)),
 ('pairing_surface_smooth',38,(240,790,785,1080)),
 ('ldos_two_inner_peaks',39,(190,110,830,525)),
 ('pairing_profiles_impurity',40,(120,278,900,476)),
 ('spectrum_scalar_a',41,(195,115,835,560)),
 ('ldos_scalar_weak',41,(270,645,755,963)),
 ('ldos_scalar_strong',42,(270,110,755,425)),
 ('ldos_scalar_middle',42,(235,575,795,935)),
 ('pairing_surface_scalar_a',43,(190,110,835,567)),
 ('spectrum_magnetic_a',43,(235,660,790,1045)),
 ('ldos_magnetic_weak',44,(270,110,755,430)),
 ('ldos_magnetic_middle',44,(230,580,795,1000)),
 ('ldos_magnetic_strong',45,(270,110,755,400)),
 ('spectrum_central_state',47,(235,115,795,507)),
 ('pairing_map_structured',48,(270,112,750,598)),
 ('pairing_surface_structured',48,(240,710,785,986)),
 ('spectrum_scalar_b',49,(235,115,795,507)),
 ('pairing_map_scalar_weak',49,(350,590,675,920)),
 ('pairing_map_scalar_strong',50,(350,112,675,440)),
 ('pairing_surface_scalar_b',50,(305,585,715,852)),
 ('spectrum_magnetic_b',51,(195,125,835,565)),
 ('pairing_surface_magnetic',51,(185,690,840,1110)),
 ('pairing_map_magnetic',52,(350,425,675,753)),
]


def main():
    target = ROOT/'diagrams'
    target.mkdir(exist_ok=True)
    records=[]
    for name,page,rectangle in REGIONS:
        source=ROOT/'reference_pages'/f'thesis_p{page:03d}.png'
        output=target/f'{name}.png'
        with Image.open(source) as original:
            original.crop(rectangle).save(output)
        records.append({'path':str(output.relative_to(ROOT)),
                        'source':str(source.relative_to(ROOT)),
                        'rectangle':list(rectangle),
                        'sha256':hashlib.sha256(output.read_bytes()).hexdigest()})
    (target/'manifest.json').write_text(json.dumps(records,indent=2)+'\n')
    print('Extracted',len(records),'diagram regions without changing plot pixels.')


if __name__ == '__main__':
    main()
