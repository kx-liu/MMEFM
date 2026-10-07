# MMEFM package development rules

MMEFM is the public R package for the Multilevel Main Effects Matrix Factor
Model.

This repository is the package repository. Research workflows, paper
reproduction code, exploratory diagnostics, historical implementations, and
large simulation outputs belong outside this repository, primarily in
`MMEFM_private`. The manuscript is maintained separately in `MMEFM_paper`.

The goal is a statistically faithful, small, readable, robust, and
CRAN-quality R package. Preserve statistical correctness before improving
implementation style, but do not preserve unnecessary implementation
complexity merely because it existed in earlier research code.

## Sources of truth

Use the manuscript in `MMEFM_paper` as the statistical specification when the
relevant method is complete and unambiguous.

Use `MMEFM_private` as a reference for established executable behaviour,
numerical conventions, regression fixtures, and research results. Do not treat
every function, helper, branch, export, or file in the private repository as
part of the desired public package.

If the manuscript and the reference implementation differ materially, do not
silently choose one. Identify the discrepancy and request or record the
intended resolution before changing statistical behaviour.

Once behaviour has been deliberately established in this public package,
focused tests and documented package contracts should protect it from
accidental change.

## Work narrowly

Start from the files and symbols directly relevant to the task.

Search for definitions and references before reading large files or unrelated
parts of the reference repositories.

Keep changes local. Do not combine a requested change with unrelated cleanup,
renaming, formatting, API redesign, or statistical modification.

Prefer one well-supported implementation over several speculative
alternatives.

Do not launch production-scale simulations automatically.

## Preserve statistical behaviour deliberately

Do not silently change statistical formulas, normalizations, identification
conventions, ranks, seeds, dimensions, numerical failure behaviour, or
returned statistical quantities.

Separate mechanical refactoring from behavioural changes.

Before a non-trivial refactor of established numerical code, identify a
deterministic regression fixture or add a focused characterization test.

Do not "improve" an estimator by adding regularization, pseudoinverses,
fallbacks, clipping, automatic rank repair, or alternative algorithms unless
that numerical policy is explicitly intended.

Robustness means detecting invalid or unsupported states clearly. It does not
mean forcing every input to produce an answer.

## Migrate selectively from MMEFM_private

Do not copy the private repository wholesale into this repository.

Move one coherent package capability at a time. Bring only the implementation
needed for the installed package and its tests.

Before migrating an internal helper, determine whether the public package
actually needs it. Inline, simplify, or omit research-era helpers when their
separate existence no longer has a clear contract.

Do not migrate paper workflows, Monte Carlo runners, benchmark scripts,
research result readers, maintenance diagnostics, historical compatibility
code, archived alternatives, or repository-specific path logic unless they are
explicitly required by the public package.

Git history and `MMEFM_private` are the archive. Do not retain obsolete code
or commented-out implementations "just in case".

## Write code a statistician would maintain

Prefer direct, ordinary, readable R code that exposes the statistical
calculation.

Do not create a helper merely to give a short piece of code a name. A clear
sequence of several lines inside its only caller is usually preferable to a
network of one-use helpers.

A helper is justified when it does at least one of the following:

- is genuinely reused;
- represents a meaningful statistical or algorithmic step;
- isolates an important numerical or validation policy;
- removes substantial complexity and has a clear independent contract.

A descriptive function name alone does not justify a helper.

Avoid wrappers that only forward arguments to another function.

Avoid helper chains whose main purpose is to hide a few elementary matrix or R
operations.

Do not introduce factories, generic utility layers, compatibility aliases,
configuration systems, extension points, or abstractions for hypothetical
future use.

Prefer cohesive functions over both giant procedures and networks of tiny
functions.

Use explicit `for` loops freely when they make a numerical algorithm easier to
read. Do not replace clear loops with nested `lapply()`, `vapply()`, `Map()`,
or similar constructs merely to make code shorter.

Avoid dense numerical one-liners.

## R style

Use `<-` for assignment.

Do not use semicolons or multiple statements on one line.

Use line breaks when they improve readability or are needed for conventional
formatting. Do not mechanically split short function calls, assignments,
conditions, or expressions across lines without a readability reason.

Use `seq_len()` and `seq_along()` instead of constructions such as `1:n` when
zero-length cases are possible.

Use `drop = FALSE` when matrix or array dimensions are part of the contract.

Do not rely on implicit vector recycling in statistical calculations.

Do not silently coerce user inputs unless the accepted coercion is obvious,
safe, and documented.

Keep mathematical names such as `Q_hat`, `J_hat`, `A1`, `Xt`, and `TT` when
they map directly to the paper. Use descriptive `snake_case` for ordinary
software concepts.

Comments should explain mathematical intent, invariants, numerical choices,
or non-obvious implementation reasons. Do not narrate obvious R syntax.

Do not use `tryCatch()` as a generic way to suppress numerical failures. Catch
only failures that are expected and intentionally handled.

Do not change the working directory, global options, environment, or random
number state unexpectedly.

Low-level deterministic numerical functions should not call `set.seed()`.
Randomness should be controlled explicitly at an appropriate public boundary
or by the caller.

## Public API

Design the public interface before exposing internal machinery.

Keep the exported API small. Export functions because they are useful to
package users, not because research scripts once called them.

Before the first CRAN release, existing interfaces inherited from research
code are provisional. Do not preserve a weak public API merely for historical
internal compatibility.

Any addition, removal, or renaming of an export must nevertheless be
deliberate and accompanied by updates to relevant callers, documentation, and
tests.

After the first public release, treat exported interfaces as
compatibility-sensitive.

Public functions should have coherent arguments and stable documented return
values. Do not expose algorithmic intermediate quantities as top-level
arguments merely because the implementation uses them.

When a fitting routine returns many related quantities, prefer a structured
classed fit object over a flat list of implementation details.

Use S3 methods when they materially improve ordinary user workflows, for
example `print()`, `summary()`, `fitted()`, or `residuals()` for model fits.

Do not create a control object or configuration layer unless the number and
nature of advanced controls genuinely justify one.

## Validation and robustness

Validate inputs at public boundaries.

Once a public boundary has established an invariant, internal functions should
normally rely on that invariant rather than repeat the same validation at every
call.

Centralize validation only for genuinely shared and non-trivial invariants. Do
not create one helper per argument or one helper per condition.

Error messages should identify the violated user-facing requirement.

Warnings should describe recoverable situations that users can reasonably act
on. Do not use warnings as routine progress messages.

Treat dimensions as explicit contracts. Where relevant, handle and test:

- zero-rank components when the model permits them;
- rank-one components;
- singleton matrix dimensions;
- R dimension dropping;
- incompatible group dimensions or ranks;
- non-finite user inputs;
- rank-deficient matrices;
- numerical boundary cases already identified by the reference implementation.

Do not replace a visible numerical failure with silent repair.

## Package boundary

Code under `R/` must justify being part of the installed package.

The installed package must not depend on `MMEFM_private`, `MMEFM_paper`, local
absolute paths, research result files, or the Git repository layout.

Paper simulations, empirical-analysis pipelines, benchmark comparisons, and
large-scale replication workflows do not belong in the installed package.

Simulation utilities may belong in the package only when they are intentionally
part of the documented user-facing functionality.

Evaluation functions used only for Monte Carlo research do not belong in the
public API.

Avoid new package dependencies when base R or an existing dependency provides a
clear and reliable implementation.

`DESCRIPTION` must contain dependencies of the installed package, not
dependencies needed only for external research workflows.

## Documentation and namespace

Use roxygen2 as the source of truth for exported API documentation and
namespace declarations.

Do not edit generated `NAMESPACE` or `.Rd` files manually unless there is a
specific reason to do so.

Document every exported function's inputs, return value, important statistical
conventions, constraints, and relevant failure modes.

Examples should be useful, deterministic where practical, and inexpensive to
run.

Keep package terminology aligned with the manuscript.

Do not document internal implementation details as if they were part of the
public contract.

## Tests

Use small deterministic tests for package behaviour.

Test statistical and dimensional invariants rather than arbitrary internal
representations.

Do not compare eigenvectors or factor coordinates naively when sign, rotation,
or equivalent representations are unidentified. Test the relevant loading
spaces, reconstructed components, objective values, or other identified
quantities instead.

Bug fixes should add a focused regression test when feasible.

Do not weaken tests or enlarge tolerances merely to make a refactor pass.

Do not put production Monte Carlo experiments in package tests.

## Validation workflow

During development, run the smallest deterministic test that covers the
changed behaviour.

Do not repeatedly run the full package check after every small edit.

After a coherent change is complete, run the relevant package tests and
regenerate documentation when necessary.

Before a release candidate, validate the package built from the source tarball,
including `R CMD check --as-cran` and the other intended release-platform
checks.

If validation reveals an unrelated pre-existing issue, report it separately
rather than expanding the current task automatically.

## Agent output

Keep implementation reports concise.

Report the files changed, the behavioural effect, the validation performed,
and any unresolved issue.

Do not restate the repository, manuscript, or these rules unless it is needed
to explain a decision.
