// to compile: first install (or just compile) htslib: https://www.htslib.org/download/
// then compile: suppose the htslib folder is in htslib-1.20
// method 1 (need to NOT change the htslib-1.20 folder position):
//    gcc -Wall -O2 -o vcf2table2 vcf2table.c -Ihtslib-1.20/htslib -Lhtslib-1.20 -lhts
// method 2 (lzma is not necssary, depending on the system configure result when compiloing htslib): 
//    gcc -O2 -Wall -g  -o vcf2table vcf2table.c htslib-1.20/libhts.a -lbz2 -lz -lm -llzma

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "htslib-1.20/htslib/hts.h"
#include "htslib-1.20/htslib/vcf.h"
#include "htslib-1.20/htslib/kstring.h"
#include "htslib-1.20/htslib/bgzf.h"

#define PACKAGE_VERSION "1.0"

int main(int argc, char *argv[]) {
    int c;
    char *mode = "wu";
    char *outfile = "-"; // stdout
    while ((c = getopt(argc, argv, "o:z")) != -1)
        switch (c) {
            case 'o':
                outfile = optarg;
                break;
            case 'z':
                mode = "w";
                break;
            case '?':
                if (optopt == 'o')
                    fprintf (stderr, "Option -%c requires an argument (output file name).\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,
                        "Unknown option character `\\x%x'.\n",
                        optopt);
                return 1;
            default:
                abort ();
      }

    if (argc - optind != 1) {
        fprintf(stderr, "\n%s version %s\n", argv[0], PACKAGE_VERSION);
        fprintf(stderr, "\nUsage: %s [options] <in.vcf(.gz)>\n\n", argv[0]);
        fprintf(stderr, "Options:\n    -o <file>     output file name [stdout]\n");
        fprintf(stderr, "\n");
        return 1;
    }

    const char *vcf_file = argv[optind];
    // FILE *output = stdout;
    // if (outfile != NULL) output = fopen(outfile, "w");
    BGZF *output = 0;
    // if (outfile != NULL) output = bgzf_open(outfile, mode);
    output = bgzf_open(outfile, mode);

    // Open the VCF file
    vcfFile *vcf = bcf_open(vcf_file, "r");
    if (!vcf) {
        fprintf(stderr, "Error opening VCF file\n");
        return 1;
    }

    // Read the VCF header
    bcf_hdr_t *hdr = bcf_hdr_read(vcf);
    if (!hdr) {
        fprintf(stderr, "Error reading VCF header\n");
        bcf_close(vcf);
        return 1;
    }

    // Print sample names (header of the table)
    kstring_t lineout = { 0, 0, NULL };
    kputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL", &lineout);
    // fprintf(output, "CHROM\tPOS\tID\tREF\tALT\tQUAL");
    int nsmpl = bcf_hdr_nsamples(hdr); // number of samples
    for (int i = 0; i < nsmpl; i++) {
        // fprintf(output, "\t%s", hdr->samples[i]);
        kputc('\t', &lineout);
        kputs(hdr->samples[i], &lineout);
    }
    // fprintf(output, "\n");
    kputc('\n', &lineout);
    ssize_t put;
    put = bgzf_write(output, lineout.s, lineout.l);
    if (put < (ssize_t) lineout.l) {
        fprintf(stderr, "Error writing : %s\n", strerror(errno));
        return -1;
    }
    lineout.l = 0;

    // Iterate over the VCF records
    bcf1_t *rec = bcf_init();
    int ngt_arr = 0, *gt_arr = NULL;
    while (bcf_read(vcf, hdr, rec) == 0) {;
        if (bcf_unpack(rec, BCF_UN_STR) != 0)
            return -1;
        // Extract basic information
        const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
        int pos = rec->pos + 1;  // VCF is 1-based
        const char *id = rec->d.id;
        const char *ref = rec->d.allele[0];
        float qual = rec->qual;

        // Print basic information
        // fprintf(output, "%s\t%d\t%s\t%s\t", chrom, pos, id, ref);
        put = ksprintf(&lineout, "%s\t%d\t%s\t%s\t", chrom, pos, id, ref);

        // Print all alternative alleles
        // fprintf(output, "%s", rec->d.allele[1]);
        kputs(rec->d.allele[1], &lineout);
        if (rec->n_allele > 2){
            for (int i = 2; i < rec->n_allele; i++) {
                // fprintf(output, ",%s", rec->d.allele[i]);
                kputc(',', &lineout);
                kputs(rec->d.allele[i], &lineout);
            }
        }
        if ( bcf_float_is_missing(qual) ) ksprintf(&lineout, "\t.");
        else ksprintf(&lineout, "\t%.1f", qual);
        // Extract alleles
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt <= 0) continue;  // No genotype information
        int max_ploidy = ngt/nsmpl;
        if (max_ploidy != 2) {
            fprintf(stderr, "your ploidy is %d\n", max_ploidy);
            perror("FIXME: currently only for diploid\n");
        }
        // Print genotypes for each sample
        for (int i = 0; i < nsmpl; i++) {
            int32_t *ptr = gt_arr + i*max_ploidy;
            if ( bcf_gt_is_missing(ptr[0]) ) ksprintf(&lineout, "\tN");
            else {
                int allele_index1 = bcf_gt_allele(ptr[0]);
                int allele_index2 = bcf_gt_allele(ptr[1]);
                if (allele_index1 == allele_index2) ksprintf(&lineout, "\t%s", rec->d.allele[allele_index1]);
                else ksprintf(&lineout, "\tH");
            }
        }
        kputc('\n', &lineout);
        put = bgzf_write(output, lineout.s, lineout.l);
        lineout.l = 0;
    }

    // Free allocated memory and close files
    free(gt_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(vcf);
    bgzf_close(output); //fclose(output);
    free(ks_release(&lineout));

    return 0;
}
