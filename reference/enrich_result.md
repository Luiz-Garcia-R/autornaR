# Print method for `enrich_result` objects

Displays a concise summary of Gene Ontology (GO) enrichment results,
including comparison information, organism, and the number of
significantly enriched terms for up- and down-regulated genes.

Returns a table of enriched Gene Ontology (GO) terms for a given
direction (up- or down-regulated genes).

## Usage

``` r
# S3 method for class 'enrich_result'
print(x, ...)

# S3 method for class 'enrich_result'
summary(object, direction = c("up", "down"), ont = NULL, ...)
```

## Arguments

- x:

  An object of class `enrich_result`.

- ...:

  Additional arguments (ignored).

- object:

  An object of class `enrich_result`.

- direction:

  Character. Which gene set to summarize: `"up"` or `"down"`.

- ont:

  Optional character specifying the ontology to return (e.g. `"BP"`,
  `"MF"`, `"CC"`). If `NULL`, returns results for all ontologies.

## Value

Invisibly returns `x`.

A data frame with enriched GO terms for the selected direction and
ontology. If `ont = NULL`, returns a list of data frames (one per
ontology).

## Details

This method is automatically called when an `enrich_result` object is
printed in the console.
