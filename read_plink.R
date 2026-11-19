


#' Read genotypes in Plink format for GWAS in RepeatABEL
#'
#' This function reads genotypes in Plink binary or Plink raw format into a
#' gwaa.data2 object that is used for genetic analysis by RepeatABEL.
#'
#' @param path Path to the files, not including file suffix. As the Plink formats
#' consist of multiple files, the convention is to refer only to the base file name
#' and suffixes like .bed, .fam, .bim or .raw will be added automatically. Note
#' that for raw format lacks map information. The function will attempt to find
#' a .map or a .bim file with the same base name to access map information. If
#' none exists, it will proceed without assigning markers to chromosomes.
#' @param format What type of Plink format to read, either "binary" (default) or
#' "raw".
#' @param ids Type of ID to use; either "iid" (default), "fid", or "concatenate".
#' Whether the animal IDs should be read from the within-family IID
#' (default), family id FID, or from the concatenation of FID and IID. The IDs
#' should correspond to the IDs used in the phenotype file. The option exists
#' because Plink format allows both within-family and family ID, and they may
#' not be unique. However, a common and recommended practice is to make the
#' within-family ID unique to each animal.
#'
#' @returns A gwaa.data2 object.
#' @export
read_plink <- function(path,
                       format = "binary",
                       ids = "iid") {
  if (format == "binary") {
    geno <- read_plink_binary(path, ids = ids)
  } else if (format == "raw") {
    geno <- read_plink_raw(path, ids = ids)
  } else {
    stop("Unknown format option.")
  }
  geno
}


read_plink_binary <- function(path,
                              ids = "iid") {

  bim <- read.table(paste0(path, ".bim"))
  colnames(bim) <- c("chr", "variant_identifier", "genetic_position", "bp", "allele1", "allele2")

  fam <- read.table(paste0(path, ".fam"))
  colnames(fam) <- c("fid", "iid", "father", "mother", "sex", "phenotype")

  ## Check ID uniqueness
  if (ids == "iid") {
    stopifnot("IDs in IID column are not unique. Try using ids = \"fid\" or ids = \"concatenate\"." =
                !any(duplicated(fam$iid)))
  } else if (ids == "fid") {
    stopifnot("IDs in FID column are not unique. Try using ids = \"iid\" or ids = \"concatenate\"." =
                !any(duplicated(fam$fid)))
  } else if (ids == "concatenate") {
  } else {
    stop("Unsupported ids option.")
  }


  geno <- read_bed(paste0(path, ".bed"),
                   n_snps = nrow(bim),
                   n_ind = nrow(fam))


  colnames(geno) <- bim$variant_identifier
  if (ids == "iid") {
    rownames(geno) <- fam$iid
  } else if (ids == "fid") {
    rownames(geno) <- fam$fid
  } else if (ids == "concatenate") {
    rownames(geno) <- paste(fam$fid, fam$iid, sep = "_")
  } else {
    stop("Unsupported ids option.")
  }

  RepeatABEL::Create_gwaa_data2(geno,
                                chromosome = as.character(bim$chr),
                                phenotypes = data.frame(id = rownames(geno)))

}


read_bed <- function(bed_file, n_snps, n_ind) {

  con <- file(bed_file, "rb")

  bed <- readBin(con, "raw",
                 n = file.info(bed_file)$size)

  close(con)

  ## Check magic number
  stopifnot(identical(bed[1:3], as.raw(c(0x6c, 0x1b, 0x01))))
  bed <- bed[-(1:3)]

  block_size <- ceiling(n_ind / 4)

  snp_block_starts <- seq(from = 1, to = n_snps * block_size, by = block_size)

  ## Check size of file
  stopifnot(length(bed) == n_snps * block_size)


  geno <- matrix(0, ncol = n_snps, nrow = n_ind)

  for (snp_ix in 1:n_snps) {

    snp <- bed[snp_block_starts[snp_ix]:(snp_block_starts[snp_ix] + block_size - 1)]

    snp_coded <- vector(mode = "list", length = length(snp))

    for (ix in 1:length(snp)) {
      b <- as.integer(snp[ix])
      snp_coded[[ix]] <- lapply(c(0, 2, 4, 6), function(s) {
        g <- bitwAnd(bitwShiftR(b, s), 3L)
        g[g == 1] <- NA
        g[g == 2] <- 1
        g[g == 3] <- 2
        g
      })
    }
    geno[, snp_ix] <- unlist(snp_coded)[1:n_ind]
  }

  geno
}


read_plink_raw <- function(path, ids = "iid") {

  raw <- read.table(paste0(path, ".raw"), header = TRUE)

  ## Check ID uniqueness
  if (ids == "iid") {
    stopifnot("IDs in IID column are not unique. Try using ids = \"fid\" or ids = \"concatenate\"." =
                !any(duplicated(raw$IID)))
  } else if (ids == "fid") {
    stopifnot("IDs in FID column are not unique. Try using ids = \"iid\" or ids = \"concatenate\"." =
                !any(duplicated(raw$FID)))
  } else if (ids == "concatenate") {
  } else {
    stop("Unsupported ids option.")
  }

  geno_matrix <- as.matrix(raw[, 7:ncol(raw)])
  colnames(geno_matrix) <- sub(colnames(geno_matrix), pattern = "_[ATGC]", replacement = "")

  if (ids == "iid") {
    rownames(geno_matrix) <- raw$IID
  } else if (ids == "fid") {
    rownames(geno_matrix) <- raw$FID
  } else if (ids == "concatenate") {
    rownames(geno_matrix) <- paste(raw$FID, raw$IID, sep = "_")
  } else {
    stop("Unsupported ids option.")
  }

  map_file <- paste0(path, ".map")
  bim_file <- paste0(path, ".bim")
  map_exists <- FALSE
  if (file.exists(map_file)) {
    message("Found map file.")
    map <- read.table(map_file)
    map_exists <- TRUE
  } else if (file.exists(bim_file)) {
    message("Found bim file.")
    map <- read.table(bim_file)
    map_exists <- TRUE
  } else {
    message("No map or bim file found. Proceeding without chromosome information.")
  }

  if (map_exists) {
    return(RepeatABEL::Create_gwaa_data2(geno_matrix,
                                         chromosome = as.character(map$V1),
                                         phenotypes = data.frame(id = rownames(geno_matrix))))
  } else {
    return(RepeatABEL::Create_gwaa_data2(geno_matrix,
                                         phenotypes = data.frame(id = rownames(geno_matrix))))
  }

}



# bed_file <- "data/simulated_missing/plink.bed"
# n_snps <- 5000
# n_ind <- 97
#
#
# bed_file <- "data/test/plink.bed"
# n_snps <- 3
# n_ind <- 6
#
#
# G <- read_bed(bed_file,
#               n_snps = n_snps,
#               n_ind = n_ind)
#
#
# gdat <- read_plink("data/test/plink", ids = "concatenate")
#
# gdat2 <- read_plink("data/simulated_missing/plink", ids = "iid")
