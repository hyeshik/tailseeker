/*
 * TailSeqCIFDumper.java
 *
 * Copyright (c) 2013 Hyeshik Chang
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

package net.sf.picard.illumina.parser;

import java.io.File;

public class TailSeqCIFDumper
{
    private static File dataDir;
    private static int laneNum;
    private static int tileNum;
    private static int numCycles;
    private static final int numChannels = 4;
    private static double signalScaleFactor;

    private static final char[] base64EncodeTable = 
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/".toCharArray();

    private static CycleIlluminaFileMap makeFileMap(final String extension) {
        final CycleIlluminaFileMap fileMap = new CycleIlluminaFileMap();
        int []cycles = new int[numCycles];

        for (int i = 0; i < numCycles; i++)
            cycles[i] = i + 1;

        fileMap.put(tileNum, new CycleFilesIterator(dataDir, laneNum, tileNum, cycles, extension));

        return fileMap;
    }

    private static CifParser makeCifParser(final String readStructure) {
        final OutputMapping outMap = new OutputMapping(new ReadStructure(readStructure));
        return new CifParser(dataDir, laneNum, makeFileMap(".cif"), outMap);
    }

    private static void run()
    {
        final String readFormat = String.format("%dT", numCycles);
        final CifParser cifParser = makeCifParser(readFormat);

        char[] intensityArray = new char[numCycles * numChannels * 2];

        int numReads;
        for (numReads = 0; cifParser.hasNext(); ++numReads) {
            final RawIntensityData rda = cifParser.next();
            int channelOffset = 0;

            for (final IntensityChannel iChan: IntensityChannel.values()) {
                short[] chanData = rda.getRawIntensities()[0].getChannel(iChan);

                for (int cycle = 0; cycle < numCycles; cycle++) {
                    short value = chanData[cycle];
                    final int arrayIndexBase = (cycle * numChannels + channelOffset) * 2;

                    value = (short)((double)value * signalScaleFactor);

                    /* Fit values in the boundary */
                    value += 255;
                    if (value >= 4096)
                        value = 4095;
                    else if (value < 0)
                        value = 0;

                    intensityArray[arrayIndexBase] = base64EncodeTable[value >> 6];
                    intensityArray[arrayIndexBase + 1] = base64EncodeTable[value & 63];
                }

                channelOffset++;
            }

            System.out.println(intensityArray);
        }
    }

    public static void main(String[] args)
    {
        laneNum = Integer.parseInt(args[1]);
        tileNum = Integer.parseInt(args[2]);
        numCycles = Integer.parseInt(args[3]);
        signalScaleFactor = Double.parseDouble(args[4]);

        dataDir = new File(args[0]);

        run();
    }
}
