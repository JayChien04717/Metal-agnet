# Project organization

Code is classified as general, qubit, resonator, or flux; layout and QComponents are separate. This structure applies to Design1, Design2, Design3, Design3_2, Design3_3, Design4, Design5, and Design6.

Shared tools are in root code/. Scripts containing design constants or cross-design dependencies stay with their Design; do not hide differences through copying or merging. Design3_2 retains its local fabrication components because they differ from shared versions.

Previously scattered files are grouped by topic under each category's legacy directory. Mixed batches and summaries remain in general/legacy. Independent runs with complete manifests retain their internal structure. results/index.csv records original paths, new paths, and content hashes; original manifests and source snapshots are not rewritten.

Removed obsolete reference scripts without active Python references and completed fixed-date batch queues; retained reusable solvers and raw results. See the [migration and cleanup record](../code/general/migration_record.json). 20261115_latest is an existing delivery snapshot and was unchanged by this reorganization.

Official parameters still require user confirmation before commit/push; the B2 candidate confirmation is retained. Moving directories does not revalidate simulations.

Validation: [migration_validation.json](../code/general/migration_validation.json). QGeometry matched before/after migration for 8 designs. Of 746 moved files, 745 had identical hashes; one Design2 manifest was updated by an active job during reorganization. Records before and after recovery analysis were retained.
