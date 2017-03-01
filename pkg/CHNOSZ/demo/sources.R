## cross-checking sources
# the reference sources
ref.source <- thermo$refs$key
# sources of elemental data
element.source <- thermo$element$source
# sources in the primary thermodynamic database
# we omit the [S92] in "HDNB78 [S92]" etc.
os1 <- gsub("\ .*", "", thermo$obigt$ref1)
os2 <- gsub("\ .*", "", thermo$obigt$ref2)
# all of the thermodynamic data sources - some of them might be NA
obigt.source <- unique(c(os1,os2))
obigt.source <- obigt.source[!is.na(obigt.source)]
# sources of protein compositions
protein.source <- thermo$protein$ref
# if the sources are all accounted for 
# these all produce character(0)
print("missing these sources for elemental properties:")
print(unique(element.source[!(element.source %in% ref.source)]))
print("missing these sources (1) for thermodynamic properties:")
print(unique(obigt.source[!(obigt.source %in% ref.source)]))
print("missing these sources for protein compositions:")
print(unique(protein.source[!(protein.source %in% ref.source)]))
# determine if all the reference sources are cited
my.source <- c(element.source,obigt.source,protein.source)
# this should produce character(0)
print("these sources are present but not cited:")
print(ref.source[!(ref.source %in% my.source)])
