# Fresh Cover20 evidence

- `independent_replay/` records the first direct local master replay of the
  original ZIP, completed in 448.59 seconds. Its additional interface report
  also records the mismatch of a separately submitted document that was not
  selected for publication. "Independent" means freshly executed locally,
  not performed by an external reviewer.
- `repository_replay/` records the release replay after lossless reconstruction
  of the split certificate in a temporary directory.
- `environment_and_inputs.json` records original source hashes and the first
  audit environment. Its `github_modified: false` describes that initial
  review, before the later explicit publication request.
- `release_checks.log` records repository input integrity and wrapper tests.
- `release_reconstruction_checks.json` records byte-for-byte restoration,
  agreement between both full replays and additional release checks.
- `SHA256SUMS` hashes this fixed evidence, excluding the checksum file itself.

The original ZIP contains no saved success reports. No old PASS flags were
used as substitutes for verification. The repository proof is the ZIP's own
Chinese text; the unrelated external Markdown is not a published proof.

From the repository root, reproduce the full chain with `make cover20-full`.
Retain new reports with a new `--report-dir` rather than overwriting these
release records. Check evidence hashes with:

```bash
cd cover20/verification
shasum -a 256 -c SHA256SUMS
```

Hashes establish byte identity, not mathematical validity or an independent
publication timestamp. Review the proof and executable acceptance together.
