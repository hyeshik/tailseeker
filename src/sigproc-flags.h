/*
 * sigproc-flags.h
 *
 * Copyright (c) 2016 Hyeshik Chang
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * - Hyeshik Chang <hyeshik@snu.ac.kr>
 */

#ifndef _SIGPROC_FLAGS_H_
#define _SIGPROC_FLAGS_H_

#define PAFLAG_POLYA_DETECTED               0x0001
#define PAFLAG_DELIMITER_HAS_MISMATCH       0x0002
#define PAFLAG_HAVE_3P_MODIFICATION         0x0004
#define PAFLAG_MEASURED_FROM_FLUORESCENCE   0x0008
#define PAFLAG_BARCODE_HAS_MISMATCHES       0x0010
#define PAFLAG_DELIMITER_IS_SHIFTED         0x0020
#define PAFLAG_DARKCYCLE_EXISTS             0x0040
#define PAFLAG_DELIMITER_NOT_FOUND          0x0080
#define PAFLAG_BALANCER_CALL_QUALITY_BAD    0x0100
#define PAFLAG_BALANCER_BIASED              0x0200
#define PAFLAG_BALANCER_SIGNAL_BAD          0x0400
#define PAFLAG_DARKCYCLE_OVER_THRESHOLD     0x0800
#define PAFLAG_CHIMERIC_TEMPLATE            0x1000  /* used by level 2 analysis */
#define PAFLAG_LIKELY_HAVE_INTACT_END       0x2000  /* used by level 2 analysis */

#endif
