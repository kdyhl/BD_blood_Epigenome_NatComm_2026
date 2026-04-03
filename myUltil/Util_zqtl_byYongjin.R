options(stringsAsFactors = FALSE)

`%c%` <- function(mat, cols) mat[, cols, drop = FALSE]
`%r%` <- function(mat, rows) mat[rows, , drop = FALSE]
`%&&%` <- function(a,b) paste(a, b, sep = '')
glue <- function(...) paste(..., sep = '')
.unlist <- function(...) unlist(..., use.names = FALSE)
.zeros <- function(n1, n2) matrix(0, n1, n2)

options(stringsAsFactors = FALSE)

.eval <- function(str) eval(parse(text = str))

write.mat <- function(mat, ...) {
    write.table(mat, col.names = FALSE, row.names = FALSE, sep = '\t',
                quote = FALSE, ...)
}

log.msg <- function(...) {
    ss = as.character(date())
    cat(sprintf('[%s] ', ss), sprintf(...), '\n', file = stderr(), sep = '')
}

.read.mat <- function(...) as.matrix(read.table(...))

read.tsv.expected <- function(in.file, col.name, col.type) {

    require(dplyr)
    require(readr)

    stopifnot(length(col.name) == length(col.type))

    ## expected columns
    expected.tab = tibble(col.name = col.name, col.type = col.type)

    ## check header
    hdr = suppressMessages(read_tsv(in.file, n_max = 1, col_names = FALSE)) %>% unlist()
    if(length(hdr) < 1) return(NULL)

    hdr.tab = tibble(col.name = hdr) %>%
        mutate(col.pos = 1:n())

    expected.hdr = hdr.tab %>% left_join(expected.tab, by = 'col.name') %>%
        mutate(col.type = if_else(is.na(col.type), '_', col.type))

    cc = expected.hdr %>% filter(col.type != '_') %>% select(col.name) %>% unlist()
    tt = expected.hdr %>% select(col.type) %>% unlist() %>% paste(collapse = '')

    if(length(cc) < nrow(expected.tab)) return(NULL)

    ret = read_tsv(in.file, col_names = cc, col_types = tt, skip = 1)

    return(ret)
}

fast.cov <- function(x, y) {
    n.obs = crossprod(!is.na(x), !is.na(y))
    ret = crossprod(replace(x, is.na(x), 0),
                    replace(y, is.na(y), 0)) / n.obs
    return(ret)
}

fast.z.cov <- function(x, y) {
    n.obs = crossprod(!is.na(x), !is.na(y))
    ret = crossprod(replace(x, is.na(x), 0),
                    replace(y, is.na(y), 0)) / sqrt(n.obs)
    return(ret)
}

## convert z-score to p-values (two-sided test)
zscore.pvalue <- function(z) {
    2 * pnorm(abs(z), lower.tail = FALSE)
}

## calculate univariate effect sizes and p-values
calc.qtl.stat <- function(xx, yy, verbose = FALSE) {

    require(dplyr)
    require(tidyr)

    .xx = scale(xx)
    .yy = scale(yy)

    ## cross-product is much faster than covariance function
    n.obs = crossprod(!is.na(.xx), !is.na(.yy))
    beta.mat = crossprod(.xx %>% rm.na.zero(), .yy %>% rm.na.zero()) / n.obs

    if(verbose) log.msg('Computed cross-products')

    ## residual standard deviation
    resid.se.mat = matrix(NA, ncol(.xx), ncol(.yy))

    for(k in 1:ncol(.yy)) {

        beta.k = beta.mat[, k]
        yy.k = .yy[, k]
        err.k = sweep(sweep(.xx, 2, beta.k, `*`), 1, yy.k, `-`)
        se.k = apply(err.k, 2, sd, na.rm = TRUE)

        if(verbose) log.msg('Residual on the column %d', k)
        resid.se.mat[, k] = se.k + 1e-8
    }

    ## organize as consolidated table
    y.cols = 1:ncol(yy)
    colnames(beta.mat) = y.cols
    colnames(n.obs) = y.cols
    colnames(resid.se.mat) = y.cols

    beta.tab = beta.mat %>%
        as.data.frame() %>%
        dplyr::mutate(x.col = 1:n()) %>%
        tidyr::gather(key = 'y.col', value = 'beta', y.cols)

    resid.se.tab = resid.se.mat %>%
        as.data.frame() %>%
        dplyr::mutate(x.col = 1:n()) %>%
        tidyr::gather(key = 'y.col', value = 'resid.se', y.cols)

    nobs.tab = n.obs %>%
        as.data.frame() %>%
        dplyr::mutate(x.col = 1:n()) %>%
        tidyr::gather(key = 'y.col', value = 'n', y.cols)

    out.tab = beta.tab %>%
        left_join(nobs.tab, by = c('x.col', 'y.col')) %>%
        left_join(resid.se.tab, by = c('x.col', 'y.col')) %>%
        dplyr::mutate(se = resid.se/sqrt(n)) %>%
        dplyr::mutate(p.val = zscore.pvalue(beta/se))

    out.tab = out.tab %>%
        mutate(x.col = as.integer(x.col)) %>%
        mutate(y.col = as.integer(y.col))

    return(out.tab)
}

fast.cor <- function(x, y) {
    x.sd = apply(x, 2, sd, na.rm = TRUE)
    y.sd = apply(y, 2, sd, na.rm = TRUE)
    ret = fast.cov(scale(x, scale = FALSE), scale(y, scale = FALSE))
    ret = sweep(sweep(ret, 1, x.sd, `/`), 2, y.sd, `/`)
    return(ret)
}

################################################################
## Negative binomial utils
adjust.size.factor <- function(xx) {
    gene.log.mean = apply(xx, 2, function(x) mean(log(x[x > 0])))
    denom = as.vector(exp(-gene.log.mean))
    size.factor = apply(t(xx), 2, function(x) median(x * denom, na.rm = TRUE))
    ret = apply(xx, 2, function(x) x / size.factor)
    return(ret)
}

nb2normal <- function(...){
    ret = log2(1/2 + ...) %>% scale() ## Voom type of transformation
    return(ret)
}

stdize.count <- function(xx) {
    .xx = xx
    .xx[xx <= 0] = NA
    xx.med = apply(.xx, 2, median, na.rm = TRUE)
    xx.med = pmax(xx.med, 1e-4)
    xx.scaled = sweep(xx, 2, xx.med, `/`)
    ret = xx.scaled * 50
    return(ret)
}

rm.na.zero <- function(xx) {
    return(replace(xx, is.na(xx), 0))
}

rm.zero <- function(xx) {
    return(replace(xx, xx == 0, NA))
}

trunc <- function(mat, lb = -4, ub = 4) {
    mat[mat > ub] = ub
    mat[mat < lb] = lb
    return(mat)
}

################################################################
## permutation
permute.rows.indep <- function(Y, rseed = 1667) {
    set.seed(rseed)
    n = nrow(Y)
    ret = apply(Y, 2,
                function(y) {
                    ret = matrix(NA, n, 1)
                    ret.pos = is.finite(y)
                    y.shuf = y[ret.pos]
                    y.shuf = sample(y.shuf)
                    ret[ret.pos, 1] = y.shuf
                    return(ret)
                })
    log.msg('Permuted breaking gene-gene correlation\n')
    return(ret)
}

################################################################
## Find most correlated (including zero values)
find.cor.idx <- function(Y1, Y0, n.ctrl, p.val.cutoff = 1) {

    colnames(Y1) = 1:ncol(Y1)
    colnames(Y0) = 1:ncol(Y0)

    require(dplyr)

    y01.stat = calc.qtl.stat(Y0, Y1) %>%
        dplyr::rename(y0 = x.col, y1 = y.col) %>%
        dplyr::filter(p.val < p.val.cutoff)

    ret = y01.stat %>% dplyr::group_by(y1) %>%
        dplyr::top_n(n = -n.ctrl, wt = p.val)

    return(ret$y0)
}

################################################################
## clean potential genetic signals
subset.plink <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

    require(dplyr)

    .error <- function(e) {
        log.msg('No QTL here!\n')
        return(NULL)
    }

    .subset <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

        chr.num = gsub(pattern = 'chr', replacement = '', chr) %>% as.integer()
        plink.cmd = sprintf('plink --bfile %s --make-bed --geno 0.05 --maf 0.05 --chr %d --from-bp %d --to-bp %d --out %s', plink.hdr, chr.num, plink.lb, plink.ub, glue(temp.dir, '/plink'))
        system(plink.cmd, ignore.stdout = TRUE, ignore.stderr = TRUE, show.output.on.console = FALSE)

        plink = zqtl::read.plink(glue(temp.dir, '/plink'))
        colnames(plink$BIM) = c('chr', 'rs', 'missing', 'snp.loc', 'plink.a1', 'plink.a2')
        colnames(plink$FAM) = c('fam', 'iid', 'father', 'mother', 'sex.code', '.pheno')
        plink$FAM = plink$FAM %>% mutate(iid = sapply(iid, gsub, pattern = 'GTEX-', replacement = ''))

        if(any(is.logical(plink$BIM$plink.a1))) {
            plink$BIM$plink.a1 = 'T'
        }

        if(any(is.logical(plink$BIM$plink.a2))) {
            plink$BIM$plink.a2 = 'T'
        }
        return(plink)
    }

    plink = tryCatch(.subset(plink.hdr, chr, plink.lb, plink.ub, temp.dir),
                     error = .error)
    return(plink)
}

## match evething to plink.gwas
match.plink <- function(plink.gwas, plink.eqtl) {

    if(is.null(plink.gwas)) return(NULL)
    if(is.null(plink.eqtl)) return(NULL)

    ret.gwas = plink.gwas
    ret.eqtl = plink.eqtl

    gwas.bim = plink.gwas$BIM %>%
        mutate(gwas.x.pos = 1:n()) %>%
        rename(gwas.plink.a1 = plink.a1,
               gwas.plink.a2 = plink.a2) %>%
        select(-missing)

    qtl.bim = plink.eqtl$BIM %>%
        mutate(qtl.x.pos = 1:n()) %>%
        rename(qtl.plink.a1 = plink.a1,
               qtl.plink.a2 = plink.a2,
               qtl.rs = rs) %>%
        select(-missing)

    bim.matched = gwas.bim %>%
        left_join(qtl.bim) %>%
        na.omit()

    if(nrow(bim.matched) < 1) return(NULL)

    bim.matched = bim.matched %>%
        dplyr::filter(((gwas.plink.a1 == qtl.plink.a1) & (gwas.plink.a2 == qtl.plink.a2)) |
                      ((gwas.plink.a2 == qtl.plink.a1) & (gwas.plink.a1 == qtl.plink.a2))) %>%
        arrange(chr, snp.loc)

    if(nrow(bim.matched) < 1) return(NULL)

    ret.gwas$BIM = ret.gwas$BIM[bim.matched$gwas.x.pos, , drop = FALSE]
    ret.gwas$BED = ret.gwas$BED[ , bim.matched$gwas.x.pos, drop = FALSE]

    ret.eqtl$BIM = ret.eqtl$BIM[bim.matched$qtl.x.pos, , drop = FALSE]
    ret.eqtl$BED = ret.eqtl$BED[ , bim.matched$qtl.x.pos, drop = FALSE]

    flip.tab = ret.gwas$BIM %>% mutate(gwas.x.pos = 1:n()) %>%
        left_join(ret.eqtl$BIM %>% mutate(qtl.x.pos = 1:n()),
                  by = c('chr', 'snp.loc'),
                  suffix = c('.gwas', '.eqtl')) %>%
        filter(plink.a1.gwas != plink.a1.eqtl)

    ret.eqtl$BIM[flip.tab$qtl.x.pos, ] = ret.gwas$BIM[flip.tab$gwas.x.pos, ]

    flip.bed = ret.eqtl$BED[, flip.tab$qtl.x.pos]
    zero.idx = flip.bed <= 0.5
    two.idx = flip.bed >= 1.5
    flip.bed[two.idx] = 0
    flip.bed[zero.idx] = 2
    ret.eqtl$BED[, flip.tab$qtl.x.pos] = flip.bed

    return(list(gwas = ret.gwas, eqtl = ret.eqtl))
}

################################################################
linearize <- function(mat, digits = 4) {
    require(dplyr)
    require(tidyr)
    colnames(mat) = 1:ncol(mat)
    ret = mat %>% as.data.frame() %>% mutate(row = 1:n()) %>%
        gather(key = col, value = val, -row) %>%
        mutate(col = as.integer(col)) %>%
        mutate(val = round(val, digits))

    return(ret)
}

melt.spike.slab <- function(effect) {
    require(dplyr)

    .lodds = linearize(effect$lodds) %>% rename(lodds = val)
    .theta = linearize(effect$theta) %>% rename(theta = val)
    .theta.sd = linearize(effect$theta.var) %>% mutate(val = sqrt(val)) %>%
        rename(theta.sd = val)

    ret = .lodds %>%
        left_join(.theta) %>%
        left_join(.theta.sd)

    return(ret)
}

melt.slab <- function(effect) {
    require(dplyr)
    ret = effect$theta %>% linearize() %>% rename(theta = val) %>%
        left_join(sqrt(effect$theta.var) %>% linearize() %>% rename(theta.sd = val))

    return(ret)
}

################################################################
sort.gwas.bed <- function(.tab) {

    if(colnames(.tab)[1] == '#chr') {
        .tab = .tab %>% rename(chr = `#chr`)
    }

    if(colnames(.tab)[1] == '#CHR') {
        .tab = .tab %>% rename(chr = `#CHR`)
    }

    ret = .tab %>%
        mutate(chr = as.character(chr)) %>%
        mutate(start = as.integer(start)) %>%
        mutate(stop = as.integer(stop)) %>%
        arrange(chr, start) %>%
        rename(`#chr` = chr)

    return(ret)
}

output.bgzip <- function(.tab, out.file, overwrite = TRUE) {

    raw.file = sub('.gz', '', out.file)

    if(file.exists(out.file) && overwrite){
        unlink(out.file)
    }

    if(file.exists(raw.file)){
        unlink(raw.file)
    }

    .tab  = sort.gwas.bed(.tab)

    readr::write_tsv(.tab, raw.file)
    log.msg('Wrote RAW file : %s', raw.file)
    Rsamtools::bgzip(raw.file, dest = out.file)
    log.msg('BGZiped RAW file : %s -> %s', raw.file, out.file)

    Rsamtools::indexTabix(out.file, 'bed')

    if(file.exists(out.file)){
        unlink(raw.file)
    }
}

