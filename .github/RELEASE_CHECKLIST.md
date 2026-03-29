# Release Checklist

## Pre-Release (1-2 weeks before)

- [ ] Create release branch: `git checkout -b release/v0.x.x`
- [ ] Update `Cargo.toml` version
- [ ] Update `CHANGELOG.md` with all changes
- [ ] Update documentation links if needed
- [ ] Run full test suite: `cargo test --lib --release`
- [ ] Run benchmarks: `cargo bench --bench alignment_benchmarks`
- [ ] Check for compiler warnings: `cargo clippy --release`
- [ ] Format code: `cargo fmt`
- [ ] Create PR for release branch review

## Release Day

- [ ] Merge release PR
- [ ] Pull latest main: `git pull origin main`
- [ ] Verify CI passes on main
- [ ] Tag release: `git tag -a v0.x.x -m "Release v0.x.x"`
- [ ] Push tag: `git push origin v0.x.x`
- [ ] Create GitHub Release from tag
  - Copy CHANGELOG.md section
  - Mark as pre-release if beta/rc
  - Set release notes

## Post-Release

- [ ] Publish to crates.io: `cargo publish`
- [ ] Verify crate published: https://crates.io/crates/omics-simd
- [ ] Update documentation site (if applicable)
- [ ] Announce on relevant channels
- [ ] Create milestone for next version
- [ ] Close related issues

## Version Numbering

- **0.1.0**: Initial release
- **0.1.x**: Bug fixes and patches
- **0.2.0**: Minor features, new matrices
- **0.3.0**: GPU acceleration
- **1.0.0**: API stability guarantee

## Semantic Versioning

- **MAJOR** (X.0.0): Breaking API changes
- **MINOR** (0.X.0): New features, backward compatible
- **PATCH** (0.0.X): Bug fixes only

## Commit Message for Release

```
release: v0.x.x

* Feature 1 description
* Feature 2 description
* Bug fix description

Breaking changes:
- None (or list if MAJOR)
```
