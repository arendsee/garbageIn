MakeGI <- function(starts, stops, scaffolds, strands=NULL, metadata=NULL, seqinfo=NULL){
  if(is.null(strands)){
    strands=rep("*", length(starts))
  } else {
    strands <- gsub('\\.', '*', strands)
  }
  g <- GenomicRanges::GRanges(
    seqnames = scaffolds,
    ranges   = IRanges::IRanges(starts, stops),
    strand   = strands,
    seqinfo  = seqinfo
  )
  if(!is.null(metadata)){
    GenomicRanges::mcols(g) <- metadata
  }
  g
}

#' Load a GFF file as a GenomicFeatures object
#'
#' @export
#' @param filename GFF filename
#' @param seqinfo_ Seqinfo object
#' @return Rmonad wrapped GenomicFeatures object
load_gene_models <- function(filename, seqinfo_=NULL){
  .load_gff(filename, tags=c("ID", "Name", "Parent"), infer_id=TRUE, seqinfo_=seqinfo_)
}

# @param tags metadata tags to keep
# @param get_naked parse untagged attributes as ID if they are they only field
# @param infer_id if ID is missing, try to create one
.load_gff <- function(file, tags, get_naked=FALSE, infer_id=FALSE, seqinfo_=NULL){

  "
  Load a GFF. There are many ways this can go wrong. Below is a summary of the
  checks and transforms applied here.

  * check GFF types (character, integer or float) 
  * unify type synonyms
  
  "

  raw_gff_ <-
    {

      "Load the raw GFF file. Raise warnings if columns have incorrect types.
      Allow comments ('#') and use '.' to indicate missing data."

      readr::read_tsv(
        file,
        col_names = c(
          "seqid",
          "source",
          "type",
          "start",
          "stop",
          "score",
          "strand",
          "phase",
          "attr"
        ),
        na        = ".",
        comment   = "#",
        col_types = "ccciidcic"
      )

    } %>_% {

      for(col in c("seqid", "type", "start", "stop")){
        if(any(is.na(.[[col]])))
          stop("GFFError: Column '", col, "' may not have missing values")
      }

  } %>>% {

    "
    Unify all type synonyms. Synonymous sets:

    gene : SO:0000704
    mRNA : messenger_RNA | messenger RNA | SO:0000234
    CDS  : coding_sequence | coding sequence | SO:0000316
    exon : SO:0000147

    The SO:XXXXXXX entries are ontology terms
    "

    gene_synonyms <- 'SO:0000704'
    mRNA_synonyms <- c('messenger_RNA', 'messenger RNA', 'SO:0000234')
    CDS_synonyms  <- c('coding_sequence', 'coding sequence', 'SO:0000316')
    exon_synonyms <- 'SO:0000147'

    .$type <- ifelse(.$type %in% gene_synonyms, 'gene', .$type)
    .$type <- ifelse(.$type %in% mRNA_synonyms, 'mRNA', .$type)
    .$type <- ifelse(.$type %in% CDS_synonyms,  'CDS',  .$type)
    .$type <- ifelse(.$type %in% exon_synonyms, 'exon', .$type)

    .

  } %>>% {

    "
    Replace transcript and coding_exon (and their synonyms) with mRNA and exon,
    respectively. This is not formally correct, but is probably the right thing
    to do. Since this is questionable, a warning is emitted if any replacements
    are made.

    The following nearly synonymous sets are merged:

    mRNA : transcript | SO:0000673
    exon : SO:0000147 | coding_exon | coding exon | SO:0000195
    "

    mRNA_near_synonyms <- c('transcript', 'SO:0000673')
    exon_near_synonyms <- c('SO:0000147', 'coding_exon', 'coding exon', 'SO:0000195')

    if(any(.$type %in% mRNA_near_synonyms)){
        .$type <- ifelse(.$type %in% mRNA_near_synonyms, 'mRNA', .$type)
        warning("Substituting transcript types for mRNA types, this is probably OK")
    }

    if(any(.$type %in% exon_near_synonyms)){
        .$type <- ifelse(.$type %in% exon_near_synonyms, 'exon', .$type)
        warning("Substituting transcript types for exon types, this is probably OK")
    }

    .

  }


  tags_ <- tags %v>%
  {

    "
    Internal. Setup and check tag list.

    1) add ID to tag list if we need to infer ID
    2) sets a temporary tag for untagged entry
    3) assert at least one tag is pressent (otherwise nothing would be done)
    "

    if(infer_id && (! "ID" %in% .)){
      . <- c("ID", .)
    }
    if((get_naked || infer_id) && (! ".U" %in% .)){
      . <- c(., ".U")
    }

    if(length(.) == 0){
      stop("No tags selected for extraction")
    }

    .

  }

  sources_ <- raw_gff_ %>>% {

    unique(as.character(.[[2]]))

  }

  raw_gff_ %>>% {

    "Extract the attribute column"

    .[[9]]

  } %>>% {

    "Split attribute column into individual fields; expressed as a dataframe
    with columns [order | ntags | tag | value]. `order` records the original
    ordering of the GFF file (which will be lost). `ntags` is a count of the
    total number of tags for a GFF row; it is a temporary column. Untagged
    values are given the temporary tag '.U', e.g. 'gene01' -> '.U=gene01'."

    tibble::data_frame(
      attr  = stringr::str_split(., ";"),
      order = seq_len(length(.))
    )                                                                                  %>>%
      dplyr::mutate(ntags = sapply(attr, length))                                      %>>%
      tidyr::unnest(attr)                                                              %>>%
      dplyr::mutate(attr = ifelse(grepl('=', attr), attr, paste(".U", attr, sep="="))) %>>%
      tidyr::separate_(
        col   = "attr",
        into  = c("tag", "value"),
        sep   = "=",
        extra = "merge"
      )

  } %>% rmonad::funnel(tags=tags_, src=sources_) %*>% {

    "Catch an AUGUSTUS shenanigan - using the tag 'Other' where 'Parent' should
    be used to link an mRNA to its gene"

    if('Other' %in% .$tag){
      if('AUGUSTUS' %in% src){
        warning("Replacing the tag 'Other' with 'Parent'. This file appears to be
        an AUGUSTUS-produced GFF, and AUGUSTUS uses the tag Other to refer to a
        Parent relationship. To avoid this warning, you may replace Other with
        Parent in the GFF files. You should confirm that Other is being used only
        to link to parents.")
        .$tag <- ifelse(.$tag == "Other", "Parent", .$tag)
      } else {
        warning("This GFF file contains the tag 'Other'. This may be fine. But
        some programs, like AUGUSTUS, use this tag to refer to a Parent. Your
        file does not appear to be from AUGUSTUS (based on source column),
        however you may want to double check the GFF. If it is using the tag
        Other as a Parent link, just replace Other with Parent and rerun fagin.
        If this file is from AUGUSTUS, you may want to change the source column
        to 'AUGUSTUS', this will allow fagin to make smarter choices in some
        circumstances.")
      }
    }

    .

  } %>% rmonad::funnel(tags=tags_) %*>% {

    "Ignore any tags other than the specified ones"

    dplyr::filter(., tag %in% tags)

  } %>_% {

    "Assert there are no commas in the extracted attribute values. These are
    legal according to the GFF spec, but I do not yet support them."

    if(any(grepl(",", .$value))){
      stop("GFFError: commas not supported in attribute tags")
    }

  } %>% rmonad::funnel(tags=tags_) %*>% {

    "Give each tag its own column"

    if(nrow(.) > 0){
      . <- tidyr::spread(., key="tag", value="value")
      for(tag in tags){
        if(!tag %in% names(.)){
          .[[tag]] <- NA_character_
        }
      }
      .
    } else {
      .$tag   = NULL
      .$value = NULL
      for(tag in tags){
        .[[tag]] <- character(0)
      }
      .
    }

  } %>>% {

    "Consider features with the parent '-' to be missing"

    if("Parent" %in% names(.)){
      .$Parent <- ifelse(.$Parent == "-", NA, .$Parent)
    }
    .

  } %>% rmonad::funnel(infer_id=infer_id) %*>% {

    "If no ID is given, but there is one untagged field, and if there are no
    other fields, then cast the untagged field as an ID. This is needed to
    accommodate the reprehensible output of AUGUSTUS."

    if(infer_id && ".U" %in% names(.)){
      # handle the excrement of AUGUSTUS
      #   * interpret .U as ID if no ID is given and if no other tags are present
      .[['ID']] <- ifelse(
        is.na(.[['ID']])            & # is this feature has no ID attribute
          (!is.na(.[['.U']]))  & # but it does have an attribute with no tag
          .[['ntags']] == 1,   # and if this untagged attribute is the only attribute
        .[['.U']],       # if so, assign the untagged attribute to ID
        .[['ID']]        # if not, just use the current ID
      )
    }

    .

  } %>% rmonad::funnel(gff=raw_gff_) %*>% {

    "Merge the attribute columns back into the GFF, remove temporary columns."

    gff$order <- seq_len(nrow(gff))

    merge(., gff, all=TRUE) %>%
      dplyr::arrange(order) %>%
      dplyr::select(-.U, -order, -ntags, -attr)

  } %>_% {

    "
    Assert the parent/child relations are correct
    "

    if(all(c("ID", "Parent") %in% names(.))){
      parents <- subset(., type %in% c("CDS", "exon"))$Parent
      parent_types <- subset(., ID %in% parents)$type

      if(any(parent_types == "gene"))
        warning("Found CDS or exon directly inheriting from a gene, this may be fine.")

      if(! all(parent_types %in% c("gene", "mRNA"))){
        offenders <- parent_types[!(parent_types %in% c("gene", "mRNA"))]
        msg <- "Found CDS or exon with unexpected parent: [%s]"
        warning(sprintf(msg, paste0(unique(offenders), collapse=", ")))
      }

      if( any(is.na(parents)) )
        stop("Found CDS or exon with no parent")

      if(! any(duplicated(.$ID, incomparables=NA)))
        warning("IDs are not unique, this is probably bad")
    }

  } %>>% {

    "Load GFF into a GenomicRanges object"

    gi <- MakeGI(
      starts    = .$start,
      stops     = .$stop,
      scaffolds = .$seqid,
      strands   = .$strand,
      metadata  = .[,c('ID','Name','Parent')]
    )

    GenomicRanges::mcols(gi)$type <- .$type

    gi

  } %>% rmonad::funnel(si=seqinfo_) %*>% {

    "Set the Seqinfo"

    # Will fail loadly if any sequence in `.` is not in the seqinfo file
    # This step is needed since seqinfo() will not create new levels.
    GenomeInfoDb::seqlevels(.) <- unique(GenomeInfoDb::seqnames(si))
    GenomicRanges::seqinfo(.) <- si 

    .

  } %>>% {

    "
    From the GenomicRanges object, create a transcript database as a TxDb object
    "

    meta <- GenomicRanges::mcols(.)

    # NOTE: This cannot just be `is_trans <- meta$type == "mRNA"` because some
    # exons are recorded as direct children of a "gene" feature. So I list as a
    # transcript anything that is the parent of an exon or CDS.
    is_trans <- meta$ID %in% meta$Parent[meta$type %in% c("exon", "CDS")]

    GenomicRanges::mcols(.)$type = ifelse(is_trans, "mRNA", meta$type)

    GenomicFeatures::makeTxDbFromGRanges(.)

  }

}
