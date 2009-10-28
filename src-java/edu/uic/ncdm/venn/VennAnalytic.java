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

import edu.uic.ncdm.venn.data.VennData;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

public class VennAnalytic {
    private int nRows;
    private int nCircles;
    private int nPolygons;
    private int nTot;
    private double stress;
    private double minStress;
    private double[] polyData;
    private double[] polyAreas;
    private double[] polyHats;
    private double[] circleData;
    private String[] circleLabels;
    private double[][] centers;
    private double[] diameters;
    private double stepsize;
    private double maxArea;
    private int totalCount;

    public VennAnalytic() {
        stepsize = .01;
        minStress = .000001;
    }

    public VennDiagram compute(VennData vd) {
        String[][] data = vd.data;
        double[] areas = vd.areas;
        boolean isAreas = vd.isAreas;

        if (isAreas)
            processAreaData(data, areas);
        else
            processElementData(data);

        computeInitialConfiguration();
        scaleDiameters(isAreas);
        scaleConfiguration();

        minimizeGlobal();

        double[][] previousCenters = new double[nCircles][2];
        copyCircles(centers, previousCenters);
        double previousStress = stress;

        minimizeLocal();
        
        if (stress > previousStress)
            copyCircles(previousCenters, centers);
        stress = computeStress();

        return collectResults();
    }

    private VennDiagram collectResults() {
        double[] colors = new double[nCircles];
        for (int j = 0; j < nCircles; j++)
            colors[j] = (double) (j + 1) / (nCircles + 1);
        double stress01 = 0;
        double stress05 = 0;
        if (nCircles > 2) {
            stress01 = Math.exp(.909 * (nCircles - 6.105)) / (1 + Math.exp(.909 * (nCircles - 6.105)));
            stress05 = Math.exp(.900 * (nCircles - 5.129)) / (1 + Math.exp(.900 * (nCircles - 5.129)));
        }
        double[] residuals = new double[nPolygons - 1];
        String[] residualLabels = new String[nPolygons - 1];
        for (int i = 1; i < nPolygons; i++) {
            residuals[i - 1] = polyAreas[i] - polyHats[i];
            char[] c = encode(i);
            String s = "";
            for (int j = 0; j < c.length; j++) {
                if (c[j] == '1')
                    s += (circleLabels[j] + "&");
            }
            s = s.substring(0, s.length() - 1);
            residualLabels[i - 1] = s;
            // System.out.println(residualLabels[i - 1] + " " + polyAreas[i] + " " + polyData[i]);
        }
        // System.out.println("stress = " + stress + ", stress01 = " + stress01 + ", stress05 = " + stress05);

        return new VennDiagram(centers, diameters, polyAreas, residuals, circleLabels, residualLabels, colors,
                stress, stress01, stress05);
    }

    private void processAreaData(String[][] data, double[] areas) {
        HashMap sets = new HashMap();
        for (int i = 0; i < data.length; i++) {
            String[] s = data[i][0].split("&");
            for (int j = 0; j < s.length; j++) {
                if (!sets.containsKey(s[j])) {
                    Double cat = new Double((double) sets.size());
                    sets.put(s[j], cat);
                }
            }
        }
        circleLabels = new String[sets.size()];
        Set keys = sets.keySet();
        Iterator it = keys.iterator();
        while (it.hasNext()) {
            String key = (String) it.next();
            int j = ((Double) sets.get(key)).intValue();
            circleLabels[j] = key;
        }
        nRows = data.length;
        nCircles = sets.size();
        nPolygons = (int) Math.pow(2, nCircles);
        polyData = new double[nPolygons];
        polyAreas = new double[nPolygons];
        polyHats = new double[nPolygons];
        circleData = new double[nCircles];
        centers = new double[nCircles][2];
        for (int i = 0; i < nRows; i++) {
            int[] subsets = new int[nCircles];
            String[] s = data[i][0].split("&");
            for (int j = 0; j < s.length; j++) {
                int jj = ((Double) sets.get(s[j])).intValue();
                subsets[jj] = 1;
            }
            int k = decode(subsets);
            polyData[k] = areas[i];
            for (int j = 0; j < nCircles; j++) {
                if (subsets[j] > 0)
                    circleData[j] += areas[i];
            }
        }
        for (int j = 0; j < nCircles; j++)
            maxArea = Math.max(maxArea, circleData[j]);
    }

    private void processElementData(String[][] data) {
        HashMap[] categories = new HashMap[2];
        categories[0] = new HashMap();
        categories[1] = new HashMap();
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < 2; j++) {
                if (!categories[j].containsKey(data[i][j])) {
                    Double cat = new Double((double) categories[j].size());
                    categories[j].put(data[i][j], cat);
                }
            }
        }
        circleLabels = new String[categories[1].size()];
        Set keys = categories[1].keySet();
        Iterator it = keys.iterator();
        while (it.hasNext()) {
            String key = (String) it.next();
            int j = ((Double) categories[1].get(key)).intValue();
            circleLabels[j] = key;
        }
        nRows = data.length;
        nCircles = categories[1].size();
        nPolygons = (int) Math.pow(2, nCircles);
        polyData = new double[nPolygons];
        polyAreas = new double[nPolygons];
        polyHats = new double[nPolygons];
        circleData = new double[nCircles];
        centers = new double[nCircles][2];
        int[][] subsets = new int[categories[0].size()][nCircles];
        for (int i = 0; i < nRows; i++) {
            int i1 = ((Double) categories[0].get(data[i][0])).intValue();
            int j1 = ((Double) categories[1].get(data[i][1])).intValue();
            subsets[i1][j1]++;
        }
        for (int i = 0; i < subsets.length; i++)
            updateCounts(subsets[i]);
    }

    protected void updateCounts(int[] counts) {
        int index = decode(counts);
        polyData[index]++;
        for (int j = 0; j < counts.length; j++) {
            if (counts[j] > 0)
                circleData[j]++;
        }
        nTot++;
    }

    private char[] encode(int index) {
        String s = Integer.toBinaryString(index);
        char[] c = s.toCharArray();
        int offset = nCircles - c.length;
        char[] result = new char[nCircles];
        for (int i = 0; i < offset; i++)
            result[i] = '0';
        System.arraycopy(c, 0, result, offset, c.length);
        return result;
    }

    private int decode(int[] subsets) {
        String b = "";
        for (int j = 0; j < subsets.length; j++) {
            if (subsets[j] > 0)
                b += '1';
            else
                b += '0';
        }
        return Integer.parseInt(b, 2);
    }

    private void calculateAreas() {
        totalCount = 0;
        int size = 200;
        byte[][][] bis = new byte[nCircles][size][size];
        double mins = Double.POSITIVE_INFINITY;
        double maxs = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < nCircles; i++) {
            double radius = diameters[i] / 2;
            mins = Math.min(centers[i][0] - radius, mins);
            mins = Math.min(centers[i][1] - radius, mins);
            maxs = Math.max(centers[i][0] + radius, maxs);
            maxs = Math.max(centers[i][1] + radius, maxs);
        }
        for (int i = 0; i < nCircles; i++) {
            double xi = (centers[i][0] - mins) / (maxs - mins);
            double yi = (centers[i][1] - mins) / (maxs - mins);
            double di = diameters[i] / (maxs - mins);
            int r = (int) (di * size / 2.);
            int r2 = r * r;
            int cx = (int) (xi * size);
            int cy = (int) (size - yi * size);
            for (int x = 0; x < size; x++) {
                for (int y = 0; y < size; y++) {
                    if ((x - cx) * (x - cx) + (y - cy) * (y - cy) < r2)
                        bis[i][x][y] = 1;
                }
            }
        }
        for (int x = 0; x < size; x++) {
            for (int y = 0; y < size; y++) {
                int[] counts = new int[nCircles];
                int count = 0;
                for (int j = 0; j < nCircles; j++) {
                    if (bis[j][x][y] == 1) {
                        counts[j]++;
                        count++;
                    }
                }
                if (count > 0)
                    updatePixels(counts);
            }
        }
        for (int i = 0; i < nPolygons; i++)
            polyAreas[i] = 100 * polyAreas[i] / totalCount;
    }

    private void updatePixels(int[] counts) {
        int index = decode(counts);
        polyAreas[index]++;
        totalCount++;
    }

    public void computeInitialConfiguration() {
        if (nCircles < 3) {
            fixedStart();
            return;
        }

        double[][] s = computeDistanceMatrix();
        double[] q = new double[nCircles];
        double[][] x = new double[nCircles][2];

        double additiveConstant = estimateAdditiveConstant(nCircles, s);

        computeScalarProducts(additiveConstant, nCircles, s, q);

        Eigen.eigenSymmetric(s, s, q);

        double[] min = new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        double[] max = new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
        for (int j = 0; j < 2; j++) {
            double rms = Math.sqrt(q[j]);
            if (Double.isNaN(rms))
                rms = 1;
            for (int i = 0; i < nCircles; i++) {
                x[i][j] = s[i][j] * rms;
                min[j] = Math.min(x[i][j], min[j]);
                max[j] = Math.max(x[i][j], max[j]);
            }
        }

        if (max[0] == min[0] || max[1] == min[1]) {
            fixedStart();
            return;
        }

        for (int i = 0; i < nCircles; i++) {
            centers[i][0] = .5 + .25 * x[i][0] / (max[0] - min[0]);
            centers[i][1] = .5 + .25 * x[i][1] / (max[1] - min[1]);
        }
    }

    private void fixedStart() {
        double theta = Math.PI / 2;
        double delta = 2 * Math.PI / nCircles;
        for (int i = 0; i < nCircles; i++) {
            centers[i][0] = .5 + Math.cos(theta);
            centers[i][1] = .5 + Math.sin(theta);
            theta -= delta;
        }
    }

    private double[][] computeDistanceMatrix() {
        double[][] s = new double[nCircles][nCircles];
        for (int i = 0; i < nPolygons; i++) {
            char[] c = encode(i);
            for (int j = 0; j < c.length; j++) {
                if (c[j] == '0')
                    continue;
                for (int k = j + 1; k < c.length; k++) {
                    if (c[k] == '0')
                        continue;
                    s[j][k] += polyData[i];
                    s[k][j] = s[j][k];
                }
            }
        }
        for (int i = 1; i < s.length; i++) {
            for (int j = 0; j < i; j++) {
                s[i][j] = -s[i][j];
                s[j][i] = s[i][j];
            }
        }
        return s;
    }

    private double estimateAdditiveConstant(int nPoints, double[][] s) {
        double additiveConstant = Double.NEGATIVE_INFINITY;
        for (int i = 1; i < nPoints; i++) {
            for (int j = 0; j < i; j++) {
                for (int k = 0; k < nPoints; k++) {
                    if (k != i && k != j) {
                        double cmin = s[i][j] - s[i][k] - s[j][k];
                        if (additiveConstant < cmin)
                            additiveConstant = cmin;
                    }
                }
            }
        }
        return additiveConstant;
    }

    private void computeScalarProducts(double additiveConstant, int nPoints, double[][] s, double[] q) {
        if (additiveConstant == Double.NEGATIVE_INFINITY)
            additiveConstant = 0.;
        double rms = 0.;
        for (int i = 1; i < nPoints; i++) {
            for (int j = 0; j < i; j++) {
                s[i][j] += additiveConstant;
                double dij = s[i][j] * s[i][j];
                rms += dij + dij;
                q[i] += dij;
                q[j] += dij;
            }
        }

        rms = rms / (nPoints * nPoints);
        double dsm;
        for (int i = 0; i < nPoints; i++) {
            for (int j = 0; j <= i; j++) {
                if (i == j)
                    dsm = 0.;
                else
                    dsm = s[i][j] * s[i][j];
                s[i][j] = ((q[i] + q[j]) / nPoints - rms - dsm) / 2.;
                s[j][i] = s[i][j];
            }
        }
    }

    private void scaleDiameters(boolean isAreas) {
        diameters = new double[nCircles];
        double area;
        for (int j = 0; j < nCircles; j++) {
            if (isAreas)
                area = circleData[j] / maxArea;
            else
                area = circleData[j] / nTot;
            diameters[j] = 2 * Math.sqrt(area / Math.PI / nCircles);
        }
    }

    private void scaleConfiguration() {
        double vc = 0;
        for (int j = 0; j < 2; j++) {
            double mc = 0;
            for (int k = 0; k < nCircles; k++)
                mc += centers[k][j];
            mc /= nCircles;
            for (int k = 0; k < nCircles; k++) {
                centers[k][j] -= mc;
                vc += centers[k][j] * centers[k][j];
            }
        }
        vc = 10. * Math.sqrt(vc / (2 * nCircles));
        if (vc > 0) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < nCircles; k++)
                    centers[k][j] /= vc;
            }
        }
    }

    private double computeStress() {
        /* regression through origin */
        calculateAreas();
        double xx = 0;
        double xy = 0;
        int n = polyData.length;
        double sst = 0;
        for (int i = 1; i < n; i++) {
            double x = polyData[i];
            double y = polyAreas[i];
            xy += x * y;
            xx += x * x;
            sst += y * y;
        }
        double slope = xy / xx;

        double sse = 0;
        for (int i = 1; i < n; i++) {
            double x = polyData[i];
            double y = polyAreas[i];
            double yhat = x * slope;
            polyHats[i] = yhat;
            sse += (y - yhat) * (y - yhat);
        }
        return sse / sst;
    }

    private void minimizeGlobal() {
        double lastStress = 1;
        for (int iter = 0; iter < 100; iter++) {
            recenter();
            stress = computeStress();
            moveGlobal();
            // System.out.println("global iteration, loss " + iter + " " + stress);
            if (stress < minStress || lastStress - stress < minStress)
                break;
            lastStress = stress;
        }
        recenter();
        stress = computeStress();
    }

    private void minimizeLocal() {
        double lastStress = 1;
        for (int iter = 0; iter < 100; iter++) {
            recenter();
            stress = computeStress();
            moveLocal();
            // System.out.println("local iteration, loss " + iter + " " + stress);
            if (stress < minStress || lastStress - stress < minStress)
                break;
            lastStress = stress;
        }
        recenter();
        stress = computeStress();
    }

    private void moveGlobal() {
        double[][] gradients = new double[nCircles][2];
        for (int i = 0; i < nPolygons; i++) {
            String s = Integer.toBinaryString(i);
            char[] c = s.toCharArray();
            int offset = nCircles - c.length;
            for (int j = 0; j < c.length; j++) {
                if (c[j] == '0')
                    continue;
                int jo = j + offset;
                for (int k = j + 1; k < c.length; k++) {
                    if (c[k] == '0')
                        continue;
                    int ko = k + offset;
                    double resid = polyAreas[i] - polyHats[i];
                    double dx = resid * stepsize * (centers[jo][0] - centers[ko][0]);
                    double dy = resid * stepsize * (centers[jo][1] - centers[ko][1]);
                    gradients[jo][0] += dx;
                    gradients[jo][1] += dy;
                    gradients[ko][0] -= dx;
                    gradients[ko][1] -= dy;
                }
            }
        }
        for (int i = 0; i < nCircles; i++) {
            centers[i][0] += gradients[i][0];
            centers[i][1] += gradients[i][1];
        }
    }

    private void moveLocal() {
        double[][] gradients = new double[nCircles][2];
        for (int i = 0; i < nCircles; i++) {
            centers[i][0] += stepsize;
            double xPlus = computeStress();
            centers[i][0] -= 2 * stepsize;
            double xMinus = computeStress();
            centers[i][0] += stepsize;
            if (xPlus < xMinus)
                gradients[i][0] = stepsize;
            else
                gradients[i][0] = -stepsize;

            centers[i][1] += stepsize;
            double yPlus = computeStress();
            centers[i][1] -= 2 * stepsize;
            double yMinus = computeStress();
            centers[i][1] += stepsize;
            if (yPlus < yMinus)
                gradients[i][1] = stepsize;
            else
                gradients[i][1] = -stepsize;
        }
        for (int i = 0; i < nCircles; i++) {
            centers[i][0] += gradients[i][0];
            centers[i][1] += gradients[i][1];
        }
    }

    private void recenter() {
        double cx = 0;
        double cy = 0;
        for (int i = 0; i < nCircles; i++) {
            cx += centers[i][0];
            cy += centers[i][1];
        }
        cx = cx / nCircles;
        cy = cy / nCircles;
        for (int i = 0; i < nCircles; i++) {
            centers[i][0] = .5 + centers[i][0] - cx;
            centers[i][1] = .5 + centers[i][1] - cy;
        }
    }

    private void copyCircles(double[][] a, double[][] b) {
        for (int i = 0; i < nCircles; i++)
            System.arraycopy(a[i], 0, b[i], 0, 2);
    }
}
