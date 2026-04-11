# GetAlphaDiversity

Calculate alpha diversity for an mmo object with flexible output modes.
Supported diversity modes:

- 'weighted' : functional Hill number (GetFunctionalHillNumber)

- 'unweighted' : Hill numbers on abundances (GetHillNumbers)

- 'faith' : Faith's phylogenetic diversity (GetFaithPD)

- 'richness' : simple feature richness (GetRichness)

## Usage

``` r
GetAlphaDiversity(
  mmo,
  q = 1,
  normalization = "None",
  mode = "richness",
  distance = "dreams",
  threshold = 0,
  filter_id = FALSE,
  id_list = NULL,
  filter_group = FALSE,
  group_list = NULL,
  output = c("sample_level", "group_average", "group_cumulative", "rarefied_sample"),
  group_col = "group",
  sample_col = "sample",
  pool_method = c("sum", "mean"),
  n_perm = 200,
  ci = 0.95,
  seed = NULL
)
```

## Arguments

- mmo:

  The mmo object containing feature data and metadata

- q:

  The Hill number order controlling abundance sensitivity (default: 1).
  Only applies to `mode = "weighted"` and `mode = "unweighted"`; ignored
  for `"richness"` and `"faith"`.

  - `q = 0` – richness: all detected features count equally regardless
    of abundance.

  - `q = 1` – Shannon-type: features weighted proportionally to their
    relative abundance.

  - `q = 2` – Simpson-type: dominant (high-abundance) features weighted
    more strongly.

- normalization:

  Abundance table to use. Options: 'None', 'Log', 'Meancentered', 'Z',
  'PA' (default: 'None'). Using `'PA'` forces presence/absence
  regardless of `mode` or `q`, effectively making every detected feature
  equally abundant before the Hill calculation.

- mode:

  The diversity metric to calculate. One of 'weighted', 'unweighted',
  'faith', 'richness' (default: 'richness'). Use `q` to control
  abundance sensitivity for 'weighted' and 'unweighted'.

- distance:

  Feature dissimilarity metric: 'dreams', 'm2ds', or 'cosine' (default:
  'dreams'). Required for `mode = "weighted"` and `mode = "faith"`;
  ignored otherwise.

- threshold:

  Numeric threshold used to define metabolite presence (default: 0)

- filter_id:

  A boolean indicating whether to filter the feature data by a specific
  list (default: FALSE)

- id_list:

  A list of feature names to filter the feature data by, if filter_id is
  TRUE (default: NULL)

- filter_group:

  A boolean indicating whether to filter the feature data by a specific
  group list (default: FALSE)

- group_list:

  A list of groups to filter the feature data by, if filter_group is
  TRUE (default: NULL)

- output:

  Output mode: 'sample_level', 'group_average', 'group_cumulative', or
  'rarefied_sample'

- group_col:

  Column in mmo\$metadata that defines groups (default: 'group')

- sample_col:

  Column in mmo\$metadata that defines sample IDs (default: 'sample')

- pool_method:

  How to pool abundances when combining samples: 'sum' or 'mean'
  (default: 'sum')

- n_perm:

  Integer; maximum number of permutations per rarefaction level
  (default: 200)

- ci:

  Numeric; confidence level (default: 0.95)

- seed:

  Optional integer seed for reproducibility (default: NULL)

## Value

For output != 'rarefied_sample': a data.frame. For output =
'rarefied_sample': a list with:

- summary: group-level rarefaction summary (mean, lwr, upr, n_perm_eff)

- raw: permutation-level values for each group and n_samples

## Details

Output modes control how samples are handled:

1.  'sample_level' : alpha per sample -\> returns sample, group, value

2.  'group_average' : mean alpha per group (summarize sample_level) -\>
    group, mean, sd, se, n, lwr, upr

3.  'group_cumulative' : pooled gamma per group (pool samples within
    group) -\> group, value

4.  'rarefied_sample' : sample-based rarefaction within group (subsample
    N samples, pool, compute) -\> group, n_samples, mean, lwr, upr,
    n_perm

NOTE: For outputs 3 and 4, pooling is performed by summing feature
intensities across samples.
