#---------------------------------------------------------------------
# Title:         IRIS - XLSX Return
# Author:        Brandon Monier
# Created:       2018-04-19 17:10:54 CDT
# Last Modified: 2018-05-16 17:16:54 CDT
#---------------------------------------------------------------------

getExcelMetadata <- function(
    wd, geo_ser, sam_out_01, sam_out_02, sam_out_03, sam_out_04,
    sam_out_05, sam_out_06, prot_out_01, pipe_out_01, md5_out,
    proc_out_01, proc_out_02, proc_out_03, raw_out_01, raw_out_02,
    raw_out_03, raw_out_04, raw_out_05, raw_out_06, pair_out_01, pair_out_02,
    pair_out_03, pair_out_04, solid_out_01, solid_out_02, solid_out_03,
    solid_out_04, solid_out_05, solid_out_06
) {

    # Copy, timestamp, and load metadata template
    file.copy("data/seq-template.xlsx", paste0("tmp/", wd))
    setwd(paste0("tmp/", wd))

    wd2 <- substring(wd, 5)
    meta <- paste0("metadata-", wd2, ".xlsx")

    file.rename(
        from = "seq-template.xlsx",
        to = meta
    )

    wb <- loadWorkbook(meta)


    # Metadata 01 - The first entry...
    meta_01a <- c(
        "# High-throughput sequencing metadata template (version 2.1).",
        "# All fields in this template must be completed.",
        paste(
            "# Templates containing example data are found in the METADATA",
            "EXAMPLES spreadsheet tabs at the foot of this page."
        ),
        paste(
            "# Field names (in blue on this page) should not be edited.",
            "Hover over cells containing field names to view field content",
            "guidelines."
        ),
        paste(
            "# Human data. If there are patient privacy concerns regarding",
            "making data fully public through GEO, please submit to NCBI's",
            "dbGaP (http://ww.ncbi.nlm.nih.gov/gap/) database. dbGaP has",
            "controlled access mechanisms and is an appropriate resource for",
            "hosting sensitive patient data."
        ),
        "",
        "SERIES",
        "# This section describes the overall experiment"
    )

    # Get number of contributers
    contrib <- geo_ser[grep("^geo_contrib_", names(geo_ser))]

    # Turn series input data into vector
    geo_ser <- unlist(geo_ser, use.names = FALSE)

    # Series title info
    ser_00_title <- c(
        "title",
        "summary",
        "overall design",
        rep("contributor", length(contrib)),
        "supplementary file",
        "SRA_center_name_code"
    )

    # The second chunk of metadata info...
    meta_02a <- c(
        "",
        "SAMPLES",
        paste(
            "# This section lists and describes each of the biological",
            "samples under investgation, as well as any protocols that are",
            "specific to individual Samples."
        ),
        paste(
            "# Additional \"processed data file\" or \"raw file\" columns",
            "may be included"
        )
    )

    # Sample `n` titles
    geo_samples <- lapply(seq_len(length(sam_out_01[[1]])), function(i) {
        paste("Sample", i)
    })
    geo_samples <- unlist(geo_samples)
    geo_samples <- c("Sample name", geo_samples)

    # Sample user titles
    geo_sam_title <- sam_out_01[1]
    geo_sam_title <- unlist(geo_sam_title, use.names = FALSE)
    geo_sam_title <- c("title", geo_sam_title)

    # Sample source name
    geo_sam_source <- sam_out_01[2]
    geo_sam_source <- unlist(geo_sam_source, use.names = FALSE)
    geo_sam_source <- c("source name", geo_sam_source)

    # Sample organism
    geo_sam_org <- sam_out_01[3]
    geo_sam_org <- unlist(geo_sam_org, use.names = FALSE)
    geo_sam_org <- c("organism", geo_sam_org)

    # Sample characteristics
    sam_out_02 <- lapply(sam_out_02, as.character)
    geo_sam_char <- lapply(seq_len(length(sam_out_02)), function(i) {
        c(
            paste("characteristic:", names(sam_out_02)[i]),
            sam_out_02[[i]]
        )
    })

    # Sample molecule
    geo_sam_mol <- sam_out_03[1]
    geo_sam_mol <- unlist(geo_sam_mol, use.names = FALSE)
    geo_sam_mol <- c("molecule", geo_sam_mol)

    # Sample description
    geo_sam_desc <- sam_out_03[2]
    geo_sam_desc <- unlist(geo_sam_desc, use.names = FALSE)
    geo_sam_desc <- c("description", geo_sam_desc)

    # Sample proc. data files (main)
    geo_sam_pdf_main <- sam_out_04
    for (i in seq_len(length(geo_sam_pdf_main))) {
        geo_sam_pdf_main[i] <- paste0(geo_sam_pdf_main[i], ".txt")
    }
    geo_sam_pdf_main <- c("processed data file", geo_sam_pdf_main)

    # Sample proc. data files (extras)
    test_05 <- sam_out_05
    for (i in seq_len(length(sam_out_05))) {
        sam_out_05[[i]] <- c("processed data file",sam_out_05[i])
    }

    # Sample raw files (all)
    test_06 <- sam_out_06
    for (i in seq_len(length(sam_out_06))) {
        sam_out_06[[i]] <- c("raw file", sam_out_06[i])
    }

    meta_03a <- c(
        "",
        "PROTOCOLS",
        paste(
            "# Any of the protocols below which are applicable to only a",
            "subset of the SAMPLES section instead."
        ),
        "growth protocol",
        "treatment protocol",
        "extract protocol",
        "library construction protocol",
        "library strategy"
    )

    geo_prot <- unlist(prot_out_01, use.names = FALSE)

    meta_04a <- c(
        "",
        "DATA PROCESSING PIPELINE",
        paste(
            "# Data processing steps include base-calling, alignment",
            "filtering, peak-calling, generation of normalized abundance",
            "measurements etc..."
        ),
        paste(
            "# For each step provide a description, as wll as software name,",
            "version, parameters, if applicable."
        ),
        "# Include additional steps, as necessary."
    )

    data_step <- pipe_out_01[grep("^data_step_", names(pipe_out_01))]
    geo_pipe <- unlist(pipe_out_01, use.names = FALSE)
    pipe_title_00 <- c(
        rep("data processing step", length(data_step)),
        "genome build",
        "processed data files format and content"
    )

    meta_05a <- c(
        "",
        paste(
            "# For each file listed in the \"processed data\" file columns",
            "of the SAMPLES section, provide additional information below."
        ),
        "PROCESSED DATA FILES"
    )

    geo_pdf_01 <- sam_out_04
    for (i in seq_len(length(geo_pdf_01))) {
        geo_pdf_01[i] <- paste0(geo_pdf_01[i], ".txt")
    }
    geo_pdf_01 <- unlist(geo_pdf_01, use.names = FALSE)
    geo_pdf_01 <- c("file name", geo_pdf_01)

    geo_pdf_02 <- c(
        "file type",
        rep("raw counts", length(geo_pdf_01) - 1)
    )

    md5_out <- unlist(md5_out, use.names = FALSE)
    md5_out <- c("file checksum", md5_out)

    geo_proc_01 <- unlist(proc_out_01, use.names = FALSE)
    geo_proc_02 <- unlist(proc_out_02, use.names = FALSE)
    geo_proc_03 <- unlist(proc_out_03, use.names = FALSE)

    meta_06a <- c(
        "",
        paste(
            "# For each file listed in the raw file columns of the SAMPLES",
            "section provide additional information below."
        ),
        "RAW FILES"
    )


    geo_raw_01 <- unlist(raw_out_01, use.names = FALSE)
    geo_raw_01 <- c("file name", geo_raw_01)

    geo_raw_02 <- unlist(raw_out_02, use.names = FALSE)
    geo_raw_02 <- c("file type", geo_raw_02)

    geo_raw_03 <- unlist(raw_out_03, use.names = FALSE)
    geo_raw_03 <- c("file checksum", geo_raw_03)

    geo_raw_04 <- unlist(raw_out_04, use.names = FALSE)
    geo_raw_04 <- c("instrument model", geo_raw_04)

    geo_raw_05 <- unlist(raw_out_05, use.names = FALSE)
    geo_raw_05 <- c("read length", geo_raw_05)

    geo_raw_06 <- unlist(raw_out_06, use.names = FALSE)
    geo_raw_06 <- c("single or paired-end", geo_raw_06)

    meta_07a <- c(
        "",
        paste(
            "# For paired-end experiments, list the 2 associated raw files,",
            "and provide average insert size and standard deviation,",
            "if known. For SOLiD experiments, list the 4 file names",
            "(include \"file name 3\" and \"file name 4\" columns)"
        ),
        "PAIRED-END EXPERIMENTS"
    )

    geo_pair_01 <- unlist(pair_out_01, use.names = FALSE)

    geo_pair_02 <- unlist(pair_out_02, use.names = FALSE)

    geo_pair_03 <- unlist(pair_out_03, use.names = FALSE)

    geo_pair_04 <- unlist(pair_out_04, use.names = FALSE)

    geo_solid_01 <- unlist(solid_out_01, use.names = FALSE)

    geo_solid_02 <- unlist(solid_out_02, use.names = FALSE)

    geo_solid_03 <- unlist(solid_out_03, use.names = FALSE)

    geo_solid_04 <- unlist(solid_out_04, use.names = FALSE)

    geo_solid_05 <- unlist(solid_out_05, use.names = FALSE)

    geo_solid_06 <- unlist(solid_out_06, use.names = FALSE)



    # Edit data...
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = meta_01a,
        startRow = 1,
        startCol = 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = ser_00_title,
        startRow = length(
            c(
                meta_01a
            )
        ) + 1,
        startCol = 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_ser,
        startRow = length(
            c(
                meta_01a
            )
        ) + 1,
        startCol = 2
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = meta_02a,
        startRow = length(
            c(
                meta_01a, ser_00_title
            )
        ) + 1,
        startCol = 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_samples,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a
            )
        ) + 1,
        startCol = 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_sam_title,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a
            )
        ) + 1,
        startCol = 2
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_sam_source,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a
            )
        ) + 1,
        startCol = 3
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_sam_org,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a
            )
        ) + 1,
        startCol = 4
    )
    for (i in seq_len(length(sam_out_02))) {
        writeData(
            wb,
            sheet = "METADATA TEMPLATE ",
            x = unlist(geo_sam_char[[i]], use.names = FALSE),
            startRow = length(
                c(
                    meta_01a, ser_00_title, meta_02a
                )
            ) + 1,
            startCol = 4 + i
        )
    }
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_sam_mol,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a
            )
        ) + 1,
        startCol = 4 + length(geo_sam_char) + 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_sam_desc,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a
            )
        ) + 1,
        startCol = 4 + length(geo_sam_char) + 2
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_sam_pdf_main,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a
            )
        ) + 1,
        startCol = 4 + length(geo_sam_char) + 3
    )

    if(!all(unlist(test_05, use.names = FALSE) == "")) {
        for (i in seq_len(length(sam_out_05))) {
            writeData(
                wb,
                sheet = "METADATA TEMPLATE ",
                x = unlist(sam_out_05[[i]], use.names = FALSE),
                startRow = length(
                    c(
                        meta_01a, ser_00_title, meta_02a
                    )
                ) + 1,
                startCol = 4 + length(sam_out_02) + 3 + i
            )
        }
    } else {
        sam_out_05 <- NULL
    }
    if(!all(unlist(test_06, use.names = FALSE) == "")) {
        for (i in seq_len(length(sam_out_06))) {
            writeData(
                wb,
                sheet = "METADATA TEMPLATE ",
                x = unlist(sam_out_06[[i]], use.names = FALSE),
                startRow = length(
                    c(
                        meta_01a, ser_00_title, meta_02a
                    )
                ) + 1,
                startCol = 4 + length(sam_out_02) + 3 + length(sam_out_05) +
                i
            )
        }
    } else {
        NULL
    }
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = meta_03a,
        startRow = length(
            c(meta_01a, ser_00_title, meta_02a, geo_samples)
        ) + 1,
        startCol = 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_prot,
        startRow = length(
            c(meta_01a, ser_00_title, meta_02a, geo_samples)
        ) + 4,
        startCol = 2
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = meta_04a,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a
            )
        ) + 1,
        startCol = 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = pipe_title_00,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a
            )
        ) + 1,
        startCol = 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_pipe,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a
            )
        ) + 1,
        startCol = 2
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = meta_05a,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00
            )
        ) + 1,
        startCol = 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_pdf_01,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a
            )
        ) + 1,
        startCol = 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_pdf_02,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a
            )
        ) + 1,
        startCol = 2
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = md5_out,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a
            )
        ) + 1,
        startCol = 3
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_proc_01,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a,
                geo_pdf_01
            )
        ) + 1,
        startCol = 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_proc_02,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a,
                geo_pdf_01
            )
        ) + 1,
        startCol = 2
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_proc_03,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a,
                geo_pdf_01
            )
        ) + 1,
        startCol = 3
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = meta_06a,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a,
                geo_pdf_01, geo_proc_01
            )
        ) + 1,
        startCol = 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_raw_01,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a,
                geo_pdf_01, geo_proc_01, meta_06a
            )
        ) + 1,
        startCol = 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_raw_02,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a,
                geo_pdf_01, geo_proc_01, meta_06a
            )
        ) + 1,
        startCol = 2
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_raw_03,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a,
                geo_pdf_01, geo_proc_01, meta_06a
            )
        ) + 1,
        startCol = 3
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_raw_04,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a,
                geo_pdf_01, geo_proc_01, meta_06a
            )
        ) + 1,
        startCol = 4
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_raw_05,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a,
                geo_pdf_01, geo_proc_01, meta_06a
            )
        ) + 1,
        startCol = 5
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_raw_06,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a,
                geo_pdf_01, geo_proc_01, meta_06a
            )
        ) + 1,
        startCol = 6
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = meta_07a,
        startRow = length(
            c(
                meta_01a, ser_00_title, meta_02a, geo_samples,
                meta_03a, meta_04a, pipe_title_00, meta_05a,
                geo_pdf_01, geo_proc_01, meta_06a, geo_raw_01
            )
        ) + 1,
        startCol = 1
    )


    if (length(geo_solid_01) != 0) {
        writeData(
            wb,
            sheet = "METADATA TEMPLATE ",
            x = c("filename 1", geo_pair_01, geo_solid_01),
            startRow = length(
                c(
                    meta_01a, ser_00_title, meta_02a, geo_samples,
                    meta_03a, meta_04a, pipe_title_00, meta_05a,
                    geo_pdf_01, geo_proc_01, meta_06a, geo_raw_01,
                    meta_07a
                )
            ) + 1,
            startCol = 1
        )
        writeData(
            wb,
            sheet = "METADATA TEMPLATE ",
            x = c("filename 2", geo_pair_02, geo_solid_02),
            startRow = length(
                c(
                    meta_01a, ser_00_title, meta_02a, geo_samples,
                    meta_03a, meta_04a, pipe_title_00, meta_05a,
                    geo_pdf_01, geo_proc_01, meta_06a, geo_raw_01,
                    meta_07a
                )
            ) + 1,
            startCol = 2
        )
        writeData(
            wb,
            sheet = "METADATA TEMPLATE ",
            x = c(
                "filename 3",
                rep("", length(geo_pair_01)),
                geo_solid_03
            ),
            startRow = length(
                c(
                    meta_01a, ser_00_title, meta_02a, geo_samples,
                    meta_03a, meta_04a, pipe_title_00, meta_05a,
                    geo_pdf_01, geo_proc_01, meta_06a, geo_raw_01,
                    meta_07a
                )
            ) + 1,
            startCol = 3
        )
        writeData(
            wb,
            sheet = "METADATA TEMPLATE ",
            x = c(
                "filename 4",
                rep("", length(geo_pair_01)),
                geo_solid_04
            ),
            startRow = length(
                c(
                    meta_01a, ser_00_title, meta_02a, geo_samples,
                    meta_03a, meta_04a, pipe_title_00, meta_05a,
                    geo_pdf_01, geo_proc_01, meta_06a, geo_raw_01,
                    meta_07a
                )
            ) + 1,
            startCol = 4
        )
        writeData(
            wb,
            sheet = "METADATA TEMPLATE ",
            x = c("average insert size", geo_pair_03, geo_solid_05),
            startRow = length(
                c(
                    meta_01a, ser_00_title, meta_02a, geo_samples,
                    meta_03a, meta_04a, pipe_title_00, meta_05a,
                    geo_pdf_01, geo_proc_01, meta_06a, geo_raw_01,
                    meta_07a
                )
            ) + 1,
            startCol = 5
        )
        writeData(
            wb,
            sheet = "METADATA TEMPLATE ",
            x = c("standard deviation", geo_pair_04, geo_solid_06),
            startRow = length(
                c(
                    meta_01a, ser_00_title, meta_02a, geo_samples,
                    meta_03a, meta_04a, pipe_title_00, meta_05a,
                    geo_pdf_01, geo_proc_01, meta_06a, geo_raw_01,
                    meta_07a
                )
            ) + 1,
            startCol = 6
        )
    } else {
        writeData(
            wb,
            sheet = "METADATA TEMPLATE ",
            x = c("filename 1", geo_pair_01),
            startRow = length(
                c(
                    meta_01a, ser_00_title, meta_02a, geo_samples,
                    meta_03a, meta_04a, pipe_title_00, meta_05a,
                    geo_pdf_01, geo_proc_01, meta_06a, geo_raw_01,
                    meta_07a
                )
            ) + 1,
            startCol = 1
        )
        writeData(
            wb,
            sheet = "METADATA TEMPLATE ",
            x = c("filename 2", geo_pair_02),
            startRow = length(
                c(
                    meta_01a, ser_00_title, meta_02a, geo_samples,
                    meta_03a, meta_04a, pipe_title_00, meta_05a,
                    geo_pdf_01, geo_proc_01, meta_06a, geo_raw_01,
                    meta_07a
                )
            ) + 1,
            startCol = 2
        )
        writeData(
            wb,
            sheet = "METADATA TEMPLATE ",
            x = c("average insert size", geo_pair_03),
            startRow = length(
                c(
                    meta_01a, ser_00_title, meta_02a, geo_samples,
                    meta_03a, meta_04a, pipe_title_00, meta_05a,
                    geo_pdf_01, geo_proc_01, meta_06a, geo_raw_01,
                    meta_07a
                )
            ) + 1,
            startCol = 3
        )
        writeData(
            wb,
            sheet = "METADATA TEMPLATE ",
            x = c("standard deviation", geo_pair_04),
            startRow = length(
                c(
                    meta_01a, ser_00_title, meta_02a, geo_samples,
                    meta_03a, meta_04a, pipe_title_00, meta_05a,
                    geo_pdf_01, geo_proc_01, meta_06a, geo_raw_01,
                    meta_07a
                )
            ) + 1,
            startCol = 4
        )
    }




    # Styles...
    bodyStyle1 <- createStyle(
        fgFill = "#CCFFCC"
    )
    bodyStyle2 <- createStyle(
        wrapText = FALSE
    )
    bodyStyle3 <- createStyle(
        fgFill = "#FFFF00",
        fontColour = "#0000FF",
        textDecoration = "bold"
    )
    bodyStyle4 <- createStyle(
        textDecoration = "bold",
        fgFill = "#CCFFCC"
    )
    bodyStyle5 <- createStyle(
        fontColour = "#FF0000",
        textDecoration = "bold"
    )
    bodyStyle6 <- createStyle(
        fontColour = "#0000FF",
        textDecoration = "bold"
    )
    bodyStyle7 <- createStyle(
        fontColour = "#0000FF"
    )
    addStyle(
        wb,
        sheet = 1,
        bodyStyle2,
        rows = 1:50,
        cols = 1:50,
        gridExpand = TRUE
    )
    addStyle(
        wb,
        sheet = 1,
        bodyStyle1,
        rows = c(
            1:6, 8,
            8 + length(ser_00_title) + 3,
            8 + length(ser_00_title) + 4,
            8 + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + 3,
            8 + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + 3,
            8 + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + 4,
            8 + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + 5,
            8 + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + length(meta_04a) +
            length(pipe_title_00) + 2,
            8 + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + length(meta_04a) +
            length(pipe_title_00) + length(meta_05a) + length(geo_pdf_01) +
            2,
            8 + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + length(meta_04a) +
            length(pipe_title_00) + length(meta_05a) + length(geo_pdf_01) +
            length(meta_06a) + length(geo_raw_01) + 2
        ),
        cols = 1:50,
        gridExpand = TRUE
    )
    addStyle(
        wb,
        sheet = 1,
        bodyStyle3,
        rows = c(
            length(meta_01a) + length(ser_00_title) + length(meta_02a) + 1,
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + length(meta_04a) +
            length(pipe_title_00) + length(meta_05a) + 1,
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + length(meta_04a) +
            length(pipe_title_00) + length(meta_05a) + length(geo_pdf_01) +
            length(meta_06a) + 1,
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + length(meta_04a) +
            length(pipe_title_00) + length(meta_05a) + length(geo_pdf_01) +
            length(meta_06a) + length(geo_raw_01) + length(meta_07a) + 1
        ),
        cols = 1:50,
        gridExpand = TRUE
    )
    addStyle(
        wb,
        sheet = 1,
        bodyStyle4,
        rows = 1:5,
        cols = 1,
        gridExpand = TRUE
    )
    addStyle(
        wb,
        sheet = 1,
        bodyStyle5,
        rows = c(
            7,
            length(meta_01a) + length(ser_00_title) + 2,
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + 2,
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + 2,
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + length(meta_04a) +
            length(pipe_title_00) + 3,
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + length(meta_04a) +
            length(pipe_title_00) + length(meta_05a) + length(geo_pdf_01) +
            3,
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + length(meta_04a) +
            length(pipe_title_00) + length(meta_05a) + length(geo_pdf_01) +
            length(meta_06a) + length(geo_raw_01) + 3
        ),
        cols = 1,
        gridExpand = TRUE
    )
    addStyle(
        wb,
        sheet = 1,
        bodyStyle6,
        rows = c(
            length(meta_01a) +
            1:length(head(ser_00_title, n = length(ser_00_title) - 2)),
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + 6,
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + 7,
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + 8,
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + length(meta_03a) + length(meta_04a) +
            1:length(pipe_title_00)
        ),
        cols = 1,
        gridExpand = TRUE
    )
    addStyle(
        wb,
        sheet = 1,
        bodyStyle7,
        rows = c(
            length(meta_01a) +
            length(head(ser_00_title, n = length(ser_00_title) - 2)) + 1,
            length(meta_01a) +
            length(head(ser_00_title, n = length(ser_00_title) - 2)) + 2,
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + 4,
            length(meta_01a) + length(ser_00_title) + length(meta_02a) +
            length(geo_samples) + 5
        ),
        cols = 1,
        gridExpand = TRUE
    )
    saveWorkbook(wb, meta, overwrite = TRUE)
}