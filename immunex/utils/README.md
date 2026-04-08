# Utils Package Policy

`immunex/utils` is a shrinking compatibility and helper layer.

Legacy product interfaces that are no longer part of the default public API
have been moved to `immunex/legacy/`.

## Allowed Contents

- generic file/path helpers
- plotting helpers
- cleanup and selection helpers
- compatibility adapters that bridge old code to the new package layout

## Modules That Need Extra Care

Several current files in this package are domain-heavy and should eventually be reviewed for migration into clearer owners, especially chain identification and topology interpretation code.

## Rule For New Code

Do not place new scientific or workflow-defining modules in `utils/` by default.

Before adding a file here, ask:

- is this generic enough to be shared by multiple subsystems?
- does it avoid owning scientific policy?
- would `analysis/`, `core/`, or `pipeline/` be a clearer owner?
