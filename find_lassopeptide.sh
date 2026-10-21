#Run inside directory antiSMASH
python3 -c '
import json, glob

for f in glob.glob("*/*.json"):
    with open(f) as fp:
        data = json.load(fp)
        for r in data.get("records", []):
            for feat in r.get("features", []):
                if "lassopeptide" in str(feat.get("qualifiers", {}).get("product", [])):
                    print(f.split("/")[0])
                    break
' | sort -u


