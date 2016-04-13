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

#define PAFLAG_HAVE_3P_MODIFICATION         1
#define PAFLAG_MEASURED_FROM_FLUORESCENCE   2
#define PAFLAG_MEASURED_USING_NAIVE_RULER   4
#define PAFLAG_DELIMITER_IS_SHIFTED         8
#define PAFLAG_BARCODE_HAS_MISMATCHES       16
#define PAFLAG_DARKCYCLE_EXISTS             32
#define PAFLAG_DELIMITER_NOT_FOUND          64
#define PAFLAG_DELIMITER_HAS_MISMATCH       128
#define PAFLAG_BALANCER_CALL_QUALITY_BAD    256
#define PAFLAG_BALANCER_BIASED              512
#define PAFLAG_BALANCER_SIGNAL_BAD          1024
#define PAFLAG_DARKCYCLE_OVER_THRESHOLD     2048
#define PAFLAG_NO_POLYA_DETECTED            4096

#endif
