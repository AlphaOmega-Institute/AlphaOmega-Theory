"""Every cited manuscript label must exist. Mathematical changes under an
unchanged label still require review of the associated check."""
import ast
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANUSCRIPTS = {
    "LG": HERE.parent / "manuscripts" / "lg" / "LG_Master_Unified.tex",
    "TAUM": HERE.parent / "manuscripts" / "taum" / "Toward_a_Universal_Model_Self_Contained.tex",
}


def _cited_labels():
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    cited = {}
    for path in sorted(HERE.rglob("test_*.py")):
        if path.name == Path(__file__).name:
            continue
        text = path.read_text()
        match = re.search(r"^LABELS\s*=\s*(\{.*?^\})", text, flags=re.S | re.M)
        if not match:
            continue
        for manuscript, labels in ast.literal_eval(match.group(1)).items():
            cited.setdefault(manuscript, {}).update({label: path.relative_to(HERE).as_posix() for label in labels})
    return cited


def test_cited_labels_exist_in_the_manuscripts():
    cited = _cited_labels()
    assert set(cited) <= set(MANUSCRIPTS) and cited
    missing = []
    for manuscript, labels in cited.items():
        source = MANUSCRIPTS[manuscript].read_text()
        defined = set(re.findall(r"\\label\{([^}]*)\}", source))
        missing += [f"{manuscript}:{label} (cited in {where})" for label, where in labels.items() if label not in defined]
    assert not missing, "labels cited by checks but absent from the manuscripts: " + ", ".join(missing)
