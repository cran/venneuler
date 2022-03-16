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

package edu.uic.ncdm.venn.data;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

public abstract class FileReader {
    private static Map separators = new LinkedHashMap();
    private static String separator = " ";
    private static int nCols;
    private static int nRows;
    private static boolean isAreas;
    private static boolean isBinaryMatrix;
    private static String[] labels;

    public static VennData getData(File fname) {
        String[][] data;
        double[] areas;
        String record;
        java.io.BufferedReader fin;
        identifySeparator(fname);
        nRows = 0;
        try {
            fin = new java.io.BufferedReader(new java.io.FileReader(fname));
            while (fin.readLine() != null)
                nRows++;
            fin.close();
            nRows--;

            data = new String[nRows][nCols];
            areas = new double[nRows];
            fin = new java.io.BufferedReader(new java.io.FileReader(fname));
            record = fin.readLine();         //header
            record = compressBlanks(record);
            String[] row = split(record);
            labels = new String[nCols];
            System.arraycopy(row, 0, labels, 0, nCols);
            for (int i = 0; i < nRows; i++) {
                record = fin.readLine();
                if (record == null)
                    break;
                record = compressBlanks(record);
                row = split(record);
                if (i == 0) {
                    if (row.length == 2)
                        isAreas = isDouble(row[1]);
                    else if (row.length > 2)
                        isBinaryMatrix = true;
                }
                if (isAreas) {
                    data[i][0] = row[0];
                    areas[i] = Double.parseDouble(row[1]);
                } else if (isBinaryMatrix) {
                    System.arraycopy(row, 0, data[i], 0, nCols);
                } else {
                    System.arraycopy(row, 0, data[i], 0, 2);
                }
            }
            fin.close();
            if (isBinaryMatrix)
                data = reformData(data);
            return new VennData(data, areas, isAreas);

        } catch (java.io.IOException ie) {
            System.err.println("I/O exception in getData");
            return null;
        }
    }

    private static String[][] reformData(String[][] data) {
        ArrayList d = new ArrayList(nRows);
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols - 1; j++) {
                data[i][j] = data[i][j].trim();
                if (data[i][j].equals("1")) {
                    String[] record = new String[2];
                    record[0] = labels[j];
                    record[1] = data[i][nCols - 1];
                    d.add(record);
                }
            }
        }
        String[][] result = new String[d.size()][2];
        for (int i = 0; i < d.size(); i++) {
            Object[] o = (Object[]) d.get(i);
            result[i][0] = (String) o[0];
            result[i][1] = (String) o[1];
        }
        return result;
    }

    private static boolean isDouble(String s) {
        try {
            Double.parseDouble(s);
        } catch (NumberFormatException ex) {
            return false;
        }
        return true;
    }

    private static void identifySeparator(File fname) {
        separators.put(",", new int[2]);
        separators.put(";", new int[2]);
        separators.put(":", new int[2]);
        separators.put("|", new int[2]);
        separators.put("\t", new int[2]);
        separators.put(" ", new int[2]);
        java.io.BufferedReader fin;
        try {
            fin = new java.io.BufferedReader(new java.io.FileReader(fname));
            String record;
            int n = 0;
            while ((record = fin.readLine()) != null) {
                record = compressBlanks(record);
                if (record.length() == 0) continue;
                countSeparators(record);
                deleteInvalidSeparators(n);
                n++;
                if (n > 100 || separators.size() == 0) break;
            }
            if (separators.size() == 0) {
                if (n == 1) {
                    separator = "";
                    nCols = 1;
                }
            } else {
                Object[] keys = separators.keySet().toArray();
                separator = (String) keys[0];
                nCols = 1 + ((int[]) separators.get(separator))[1];
            }
            fin.close();
        } catch (java.io.IOException ie) {
            System.err.println("I/O exception in computeVariableTypes");
        }
    }

    private static void countSeparators(String record) {
        boolean isUnquoted = true;
        for (int i = 0; i < record.length(); i++) {
            if (record.charAt(i) == '"')
                isUnquoted = !isUnquoted;
            if (isUnquoted) {
                Object[] keys = separators.keySet().toArray();
                for (int j = 0; j < keys.length; j++) {
                    String key = (String) keys[j];
                    if (record.substring(i, i + 1).equals(key)) {
                        int[] counts = (int[]) separators.get(key);
                        counts[0]++;
                        break;
                    }
                }
            }
        }
    }

    private static void deleteInvalidSeparators(int n) {
        Object[] keys = separators.keySet().toArray();
        for (int j = 0; j < keys.length; j++) {
            String key = (String) keys[j];
            int[] counts = (int[]) separators.get(key);
            if (counts[0] == 0 || n > 0 && counts[0] != counts[1]) {
                separators.remove(key);
            } else {
                counts[1] = counts[0];
                counts[0] = 0;
                separators.put(key, counts);
            }
        }
    }

    private static String[] split(String record) {
        String[] row = new String[nCols];
        boolean isUnquoted = true;
        int i0 = 0;
        int col = 0;
        for (int i = 0; i < record.length(); i++) {
            if (record.charAt(i) == '"')
                isUnquoted = !isUnquoted;
            if (isUnquoted) {
                if (record.substring(i, i + 1).equals(separator)) {
                    row[col] = record.substring(i0, i);
                    row[col] = trimSpaces(row[col]);
                    row[col] = trimQuotes(row[col]);
                    col++;
                    i0 = i + 1;
                } else if (i == record.length() - 1) {
                    row[col] = record.substring(i0, i + 1);
                    row[col] = trimSpaces(row[col]);
                    row[col] = trimQuotes(row[col]);
                }
            }
        }
        if (row[nCols - 1] == null)
            row[nCols - 1] = "";
        return row;
    }

    private static String trimSpaces(String s) {
        return s.trim();
    }

    private static String trimQuotes(String s) {
        if (s.startsWith("\""))
            s = s.substring(1, s.length());
        if (s.endsWith("\""))
            s = s.substring(0, s.length() - 1);
        return s;
    }

    private static String compressBlanks(String record) {
        return record.replaceAll(" {2,}", " ").trim();
    }

}