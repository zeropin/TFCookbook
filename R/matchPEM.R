#' @importClassesFrom GenomicRanges GenomicRanges

setGeneric("matchPEM",
           function(PEM, subject, ...) standardGeneric("matchPEM"))

#' find binding sites within some GenomicRanges or DNA sequences that match some PEM
#' 
#' @describeIn matchPEM matrix/GenomicRanges
#' @export
#' @param PEM a PEM matrix
#' @param subject a GenomicRanges
#' @param genome the genome name, e.., hg38, mm10, etc
#' @param E.cutoff the threshold to select hit
#' @example matchPEM(PEM = some.PEM, subject = some.GRanges, genome = "hg38", E.cutoff = -3)
setMethod("matchPEM", signature(PEM     = "matrix",
                                subject = "GenomicRanges"),
          function(PEM,
                   subject,
                   genome = GenomeInfoDb::genome(subject),
                   E.cutoff = -3) {
            
            PEM <- .validate_PEM_input(PEM)
            GenomicRanges::strand(subject) <- "+"
            genome <- .validate_genome_input(genome)
            seqs <- getSeq(x = genome, names = subject, as.character = TRUE)
            
            w <- 7
            
            tmp_out <- .get_motif_positions(list(-PEM), seqs, nuc_freqs = rep(0.25, 4), -E.cutoff, w)
            
            m_ix <- which(tmp_out$motif_ix == 0)  
            out <- GRanges(seqnames = seqnames(subject)[tmp_out$seq_ix[m_ix] +  1],
                           ranges = IRanges(start = start(subject[tmp_out$seq_ix[m_ix] + 1]) + tmp_out$pos[m_ix],
                                            width = ncol(PEM)),
                           strand = tmp_out$strand[m_ix]) 
              
            if(is.null(names(subject)))
              out$fragment.ID <- tmp_out$seq_ix[m_ix] +1
            else
              out$fragment.ID <- names(subject)[tmp_out$seq_ix[m_ix] +1]
              
            out$Sequence <- getSeq(x = genome, names = out, as.character = TRUE)
            out$predicted.Energy <- -tmp_out$score[m_ix]
            out <- unique(out)
            
            return(out)
          })

#' @describeIn matchPEM matrix/character
#' @export
#' @param PEM a PEM matrix
#' @param subject a list of DNA sequences
#' @param E.cutoff the threshold to select hit
#' @example matchPEM(PEM = somePEM, subject = sequences, E.cutoff = -3)
setMethod("matchPEM", signature(PEM     = "matrix",
                                subject = "character"),
          function(PEM,
                   subject,
                   E.cutoff = -3) {
            PEM <- .validate_PEM_input(PEM)
            w    <- 7
            
            if(is.null(names(subject)))
              names(subject) <- 1:length(subject)
              
            tmp_out <- .get_motif_positions(list(-PEM), subject, nuc_freqs = rep(0.25, 4), -E.cutoff, w)
              
            m_ix <- which(tmp_out$motif_ix == 0)  
            out <- GRanges(seqnames = names(subject)[tmp_out$seq_ix[m_ix] +  1],
                           ranges   = IRanges(start = tmp_out$pos[m_ix] + 1,
                                              width = ncol(PEM)),
                           strand   = tmp_out$strand[m_ix])
            
            out$Sequence <- getSeq(x = DNAStringSet(subject), names = out)
            out$predicted.Energy <- -tmp_out$score[m_ix] 
            out <- unique(out)
            
            return(as_tibble(out))
          })


#' @describeIn matchPEM matrix/DNAStringSet
#' @export
setMethod("matchPEM", signature(PEM     = "matrix",
                                subject = "DNAStringSet"),
          function(PEM,
                   subject,
                   E.cutoff = -3) {
            return(matchPEM(PEM, as.character(subject), E.cutoff))
          })




setGeneric("matchSite",
           function(site, subject, ...) standardGeneric("matchSite"))

#' @describeIn matchSite character/GenomicRanges
#' @export
#' @example
#' matchSite("GAAGCG", GRanges, mismatch = 1)
setMethod("matchSite", signature(site   = "character",
                                 subject= "GenomicRanges"),
          function(site,
                   subject,
                   genome,
                   max.mismatch = 0){
            
              PEM<-addAnchorMatrix(PEM = NULL, anchor = site, height = 1)
              out <- matchPEM(PEM,
                              subject = subject,
                              genome  = genome,
                              E.cutoff = -nchar(site)+max.mismatch+0.1)
    
              out$predicted.Energy = NULL
              return(out)
})

#' @describeIn matchSite character/DNAStringSet
#' @export
#' @example
#' matchSite("GAAGCG", some DNAStringSet, mismatch = 1)
setMethod("matchSite", signature(site   = "character",
                                 subject= "DNAStringSet"),
          function(site,
                   subject,
                   max.mismatch = 0){
            
              PEM<-addAnchorMatrix(PEM = NULL, anchor = site, height = 1)
              out <- matchPEM(PEM,
                              subject = subject,
                              E.cutoff = -nchar(site)+max.mismatch+0.1)
          
              out$predicted.Energy = NULL
              return(out)
          })


#' @describeIn matchSite character/character
#' @export
#' @example
#' matchSite("GAAGCG", "AGAAGCGTTACG", mismatch = 1)
setMethod("matchSite", signature(site   = "character",
                                 subject= "character"),
          function(site,
                   subject,
                   max.mismatch = 0){
            
              PEM<-addAnchorMatrix(PEM = NULL, anchor = site, height = 1)
              out<-matchPEM(PEM,
                            subject = subject,
                            E.cutoff = -nchar(site)+max.mismatch+0.1)
          
              out$predicted.Energy = NULL
              return(out)
          })

