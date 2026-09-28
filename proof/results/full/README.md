# The certificate

[`search.cert`](search.cert) is the complete certificate for the whole search
domain, produced by `rid search` with [`search.json`](search.json) and checked
by `rid check` with [`check.json`](check.json), which reports it complete
([`check-summary.json`](check-summary.json)).

| | |
| --- | --- |
| Records | 192,696 |
| Evaluated boxes | 385,391 |
| Deepest record | depth 49 (mean 29.6) |
| Size | 18,195,036 bytes |
| SHA-256 | `df1f8dbab68b561f22fa9daeca48153d71bde7b6fba3345d9f6b7b8077df4b5f` |
| Search | 1 min 9 s on 16 threads (about 18 processor minutes) |
| Check | 36 s on 16 threads |

Records by component: Global 157,458, Local 21,792, Domain 12,650,
Exotic 796. The same file is produced for every number of threads, and
after an interrupted search is resumed.

The records were also re-verified by the independent checker in
[`../../../audit`](../../../audit).
