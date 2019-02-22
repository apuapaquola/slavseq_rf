#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"

#include "sw.h"

MODULE = Gapafim		PACKAGE = Gapafim		

void
sw(first, second)
       char* first;
       char* second;
    PPCODE:
       sw(first, second, &sp);
