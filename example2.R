#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")

parser$add_argument("-q", "--quietly", action="store_false",
                    dest="verbose", help="Print little output")

parser$add_argument("-r", "--repeats", type="integer", default=1,
                    help="Number of repeats used for cross-validation [default %(default)s]",
                    metavar="number")

parser$add_argument("-b", "--best", action="store_false", default=FALSE,
                    help="Use best predictors found in linear regression [default]")

parser$add_argument("input_file", nargs=1, help="File to be displayed")


args <- parser$parse_args()
input_file <- args$input_file

if( file.access(input_file) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", input_file))
} else {
  input.data <- get(load(input_file))
}

print(head(input.data))
