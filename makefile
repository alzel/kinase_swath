O = "./R/objects"

all: data

data: ./output/load.Rout ./output/clean.Rout

./output/load.Rout: ./R/boot.R ./R/functions.R ./R/load.R ./data/*
	R CMD BATCH --no-save --no-restore ./R/load.R ./output/load.Routi

./output/clean.Rout: ./R/clean.R ./R/objects/*_load_.RData
	R CMD BATCH --no-save --no-restore ./R/clean.R ./output/clean.Rout

#data: ./R/*_load_.RData ./R/*_clean_.RData

#loading data
#./R/*_load_.RData: ./R/boot.R ./R/load.R ./data/*
#	R CMD BATCH --no-save --no-restore ./R/load.R ./output/load.Rout

#cleaning data
#./R/*_clean_.RData: ./R/boot.R ./R/clean.R
#	R CMD BATCH --no-save --no-restore ./R/clean.R ./output/clean.Rout



clean:
	rm ./output/*
