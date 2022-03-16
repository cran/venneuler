/*
 * VennEuler -- A Venn and Euler Diagram program.
 *
 * Copyright 2009 by Leland Wilkinson.
 *
 * The contents of this file are subject to the Mozilla Public License Version 1.1 (the "License")
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
 * for the specific language governing rights and limitations under the License.
 */

package edu.uic.ncdm.venn;

import edu.uic.ncdm.venn.data.FilePicker;
import edu.uic.ncdm.venn.data.FileReader;
import edu.uic.ncdm.venn.data.VennData;
import edu.uic.ncdm.venn.display.VennFrame;

import java.io.File;

public class Venn {

    private Venn() {}

    public static void main(String[] argv) {
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                process();
            }
        });
    }

    private static void process() {
        String fileName = FilePicker.loadFile();
        File file = new File(fileName);
        VennData dv = FileReader.getData(file);
        VennAnalytic va = new VennAnalytic();
        VennDiagram vd = va.compute(dv);
        new VennFrame(vd);
    }
}
