      PROGRAM ANNOFETCH
      CHARACTER*200 VCFFILE,OUTFILE,CHRFILE(25),BL80,LINE
      CHARACTER*20 CH20
      CHARACTER*1 BLNK,CHROM(250000000),HTAB
      INTEGER LENCHR(25)
      DATA BL80/' '/
      DATA BLNK/' '/
      CALL GETARG(1,VCFFILE)
      CALL GETARG(2,CH20)
      CALL GETARG(3,OUTFILE)
      READ(CH20,*) NDELT
C
      HTAB=CHAR(9)
      OPEN(UNIT=1,FILE='/bioinfoC/lukebria/chrom/chr-length.txt',
     *FORM='FORMATTED')
      REWIND(1)
      READ(1,'(I10,1X,A80)')(LENCHR(I),CHRFILE(I),I=1,25)
      CLOSE(1)
C
      ICHROLD=0
      OPEN(UNIT=1,FILE=VCFFILE,FORM='FORMATTED')
      REWIND(1)
      OPEN(UNIT=2,FILE=OUTFILE,FORM='FORMATTED')
      REWIND(2)
    8 READ(1,'(A)',ERR=36,END=36) LINE
      IF(LINE(1:1).EQ.'#') THEN
        WRITE(2,'(A)') 'ERR'
        GOTO 8
      ENDIF
      IJ=INDEX(LINE,HTAB)
      IF(IJ.LT.2.OR.IJ.GT.6) THEN
        WRITE(2,'(A)') 'ERR'
        GOTO 8
      ENDIF
      IF(LINE(1:1).EQ.'C'.OR.LINE(1:1).EQ.'c') THEN
        IF(IJ.EQ.5) THEN
          IF(LINE(4:4).EQ.'X'.OR.LINE(4:4).EQ.'x') THEN
            ICHRNEW=23
          ELSEIF(LINE(4:4).EQ.'Y'.OR.LINE(4:4).EQ.'y') THEN
            ICHRNEW=24
          ELSEIF(LINE(4:4).EQ.'M'.OR.LINE(4:4).EQ.'m') THEN
            ICHRNEW=25
          ELSE
            READ(LINE,'(3X,I1)',ERR=14,END=14) ICHRNEW
          ENDIF
        ELSE
          READ(LINE,'(3X,I2)',ERR=14,END=14) ICHRNEW
        ENDIF
      ELSE
        IF(LINE(1:1).EQ.'X'.OR.LINE(1:1).EQ.'x') THEN
          ICHRNEW=23
        ELSEIF(LINE(1:1).EQ.'Y'.OR.LINE(1:1).EQ.'y') THEN
          ICHRNEW=24
        ELSEIF(LINE(1:1).EQ.'M'.OR.LINE(1:1).EQ.'m') THEN
          ICHRNEW=25
        ELSE
          READ(LINE,*,ERR=14,END=14) ICHRNEW
        ENDIF
      ENDIF
      LINE(1:IJ)=BL80(1:IJ)
      READ(LINE,*,ERR=14,END=14) MIDP
      ISTART=MIDP-NDELT
      ISTOP=MIDP+NDELT
C
      IF(ICHRNEW.NE.ICHROLD) THEN
        OPEN(UNIT=11,FILE=CHRFILE(ICHRNEW),FORM='UNFORMATTED')
        REWIND(11)
        READ(11)(CHROM(I),I=1,LENCHR(ICHRNEW))
        CLOSE(11)
        ICHROLD=ICHRNEW
      ENDIF
C
      WRITE(2,'(10000A1)') (CHROM(J),J=ISTART,ISTOP)
      GOTO 8
C
   14 WRITE(2,'(A)') 'ERR'
      GOTO 8
C
   36 CLOSE(1)
      CLOSE(2)
      STOP
      END
