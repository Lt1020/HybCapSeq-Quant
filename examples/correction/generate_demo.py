"""Generate explicitly simulated paired concentrations with analytical truth."""

import csv
from pathlib import Path


def main():
    rows = []
    for parent in range(4):
        for extraction in range(2):
            for dilution in range(3):
                for method, coefficient in [
                    ("depth_calibration", 10),
                    ("sequencing_yield", 5),
                ]:
                    truth = (parent + 1) * (extraction + 1) * 100
                    rows.append(
                        dict(
                            batch="SIMULATED_BATCH",
                            sample_id=f"SIM_P{parent}_E{extraction}_D{dilution}",
                            parent=f"SIM_P{parent}",
                            extraction_id=f"SIM_P{parent}_E{extraction}",
                            method=method,
                            target="SIMULATED_TARGET",
                            seq_extract_copies_ul=truth / coefficient,
                            dpcr_extract_copies_ul=truth,
                            complete_in_range=True,
                            dpcr_comparator_eligible=True,
                            model_accepted=True,
                            data_origin="simulated_analytical_truth",
                        )
                    )
    path = Path(__file__).parent / "simulated_pairs.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    print(path)


if __name__ == "__main__":
    main()
