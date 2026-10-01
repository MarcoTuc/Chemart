;redcode-94
;name Conformance
;author Chemart
        ORG    start
start   MOV.X  xsrc, xdst
        SEQ.I  xdst, xwant
        JMP    fail
        MOV.I  isrc, idst
        SEQ.I  idst, iwant3
        JMP    fail
        SNE.I  nsrc, ndst
        JMP    fail
        MOV.AB #9, abdst
        SEQ.I  abdst, abwant
        JMP    fail
        MOV.AB #9, >iptr
        SEQ.I  iptr, iwant
        JMP    fail
        SEQ.I  itgt, iwant2
        JMP    fail
        MOV.AB #7, <dptr
        SEQ.I  dptr, dwant
        JMP    fail
        SEQ.I  dtgt, dwant2
        JMP    fail
        ADD.AB #1, wrap
        SEQ.I  wrap, wwant
        JMP    fail
        DJN.B  fail, dcount
        SEQ.I  dcount, dcwant
        JMP    fail
        SLT.AB #1, sltv
        JMP    fail
        SLT.AB #9, sltv
        JMP    done
        JMP    fail
done    JMP    done
fail    DAT.F  #0, #0
xsrc    DAT.F  #1, #2
xdst    DAT.F  #0, #0
xwant   DAT.F  #2, #1
isrc    SPL.B  #7, <3
idst    DAT.F  #0, #0
iwant3  SPL.B  #7, <3
nsrc    DAT.F  #0, #1
ndst    DAT.F  $0, #1
abdst   DAT.F  #5, #6
abwant  DAT.F  #5, #9
iptr    DAT.F  #0, $2
iwant   DAT.F  #0, $3
itgt    DAT.F  #0, #0
iwant2  DAT.F  #0, #9
dptr    DAT.F  #0, $3
dwant   DAT.F  #0, $2
dtgt    DAT.F  #0, #0
dwant2  DAT.F  #0, #7
wrap    DAT.F  #0, #799
wwant   DAT.F  #0, #0
dcount  DAT.F  #0, #1
dcwant  DAT.F  #0, #0
sltv    DAT.F  #0, #5
        END
