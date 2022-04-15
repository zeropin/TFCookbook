#' @importClassesFrom GenomicRanges GenomicRanges

setGeneric("matchPEM",
           function(PEM, subject, ...) standardGeneric("matchPEM"))

#' @describeIn matchPEM matrix/GenomicRanges
#' @export
setMethod("matchPEM", signature(PEM     = "matrix",
                                subject = "GenomicRanges"),
          function(PEM,
                   subject,
                   genome = GenomeInfoDb::genome(subject),
                   out = c("positions", "scores"),
                   E.cutoff = -3) {
            out <- match.arg(out)
            PEM <- .validate_PEM_input(PEM)
            GenomicRanges::strand(subject) <- "+"
            genome <- .validate_genome_input(genome)
            seqs <- getSeq(x = genome, names = subject, as.character = TRUE)
            
            w <- 7
              
            
            if (out == "positions") {
              tmp_out <- .get_motif_positions(list(-PEM), seqs, nuc_freqs = rep(0.25, 4), -E.cutoff, w)
            
              m_ix <- which(tmp_out$motif_ix == 0)  
              out <- GRanges(seqnames = seqnames(subject)[tmp_out$seq_ix[m_ix] +  1],
                           ranges = IRanges(start = start(subject[tmp_out$seq_ix[m_ix] + 1]) + tmp_out$pos[m_ix],
                                            width = ncol(PEM)),
                           strand = tmp_out$strand[m_ix]) 
              
              names(out) <- names(PEM)
              out$Sequence <- getSeq(x = genome, names = out, as.character = TRUE)
              out$predicted.Energy <- -tmp_out$score[m_ix]
              out <- unique(out)
            }
            else if(out == "scores"){
              tmp_out <- .get_motif_positions(list(-PEM), seqs, nuc_freqs = rep(0.25, 4), -E.cutoff, w)
              
              m_ix <- which(tmp_out$motif_ix == 0)
              
              out  <- -tmp_out$score[m_ix]
            }
            
            return(out)
          })

#' @describeIn matchPEM matrix/DNAStringSet
#' @export
setMethod("matchPEM", signature(PEM     = "matrix",
                                subject = "DNAStringSet"),
          function(PEM,
                   subject,
                   out = c("positions", "scores"),
                   E.cutoff = -3) {
            out  <- match.arg(out)
            seqs <- as.character(subject)
            w    <- 7
            
            if(is.null(names(subject)))
              names(subject) <- 1:length(subject)
            
            if (out == "positions") {
              tmp_out <- .get_motif_positions(list(-PEM), seqs, nuc_freqs = rep(0.25, 4), -E.cutoff, w)
              
              m_ix <- which(tmp_out$motif_ix == 0)  
              out <- GRanges(seqnames = names(subject)[tmp_out$seq_ix[m_ix] +  1],
                             ranges   = IRanges(start = tmp_out$pos[m_ix] + 1,
                                                width = ncol(PEM)),
                             strand   = tmp_out$strand[m_ix]) 
              
              out$Sequence <- as.character(getSeq(subject, out))
              out$predicted.Energy <- -tmp_out$score[m_ix]
              out <- unique(out)
            }
            else if(out == "scores"){
              tmp_out <- .get_motif_positions(list(-PEM), seqs, nuc_freqs = rep(0.25, 4), -E.cutoff, w)
              
              m_ix <- which(tmp_out$motif_ix == 0)
              
              out  <- -tmp_out$score[m_ix]
            }
            
            return(out)
          })


#' @describeIn matchPEM matrix/character
#' @export
setMethod("matchPEM", signature(PEM     = "matrix",
                                subject = "character"),
          function(PEM,
                   subject,
                   out = c("positions", "scores"),
                   E.cutoff = -3) {
            return(matchPEM(PEM, DNAStringSet(subject), out, E.cutoff))
          })
