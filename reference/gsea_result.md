# Print method for gsea_result objects

Displays a formatted summary of a `gsea_result` object, including
contrast, organism, and top enriched pathways.

Returns a table of top enriched pathways from a `gsea_result` object.

## Usage

``` r
# S3 method for class 'gsea_result'
print(x, ...)

# S3 method for class 'gsea_result'
summary(object, ...)
```

## Arguments

- x:

  An object of class `gsea_result`.

- ...:

  Further arguments (ignored).

- object:

  An object of class `gsea_result`.

## Value

Invisibly returns `x`.

A data.frame containing the top GSEA results.

## Details

gsea_result object

An object returned by
[`rna_gsea()`](https://luiz-garcia-r.github.io/autornaR/reference/rna_gsea.md)
containing GSEA results, rankings, and plots.

## Slots

- ranking:

  Named vector

- gsea_full:

  Full table

- gsea_top:

  Filtered table
