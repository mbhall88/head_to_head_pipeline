import sys


def main():
    truth = sys.argv[1]
    classification = sys.argv[2]

    truth_map = dict()
    with open(truth) as stream:
        for line in map(str.rstrip, stream):
            if not line:
                continue
            lineage, sample = line.split("\t")
            truth_map[sample] = lineage.replace("L", "")

    with open(classification) as csvfile:
        _ = next(csvfile)
        for row in map(str.rstrip, csvfile):
            if not row:
                continue
            fields = row.split(",")
            sample = fields[0]
            major_lineage = fields[1]
            if major_lineage == truth_map[sample]:
                result = "PASS"
            else:
                result = "FAIL"
            print(f"{sample}: [{result}]")
            if result == "FAIL":
                print(
                    f"Got {major_lineage}, expected {truth_map[sample]}",
                    file=sys.stderr,
                )


if __name__ == "__main__":
    main()
