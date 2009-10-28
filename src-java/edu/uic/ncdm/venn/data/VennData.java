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
 *
 * CHANGES:
 *   Simon Urbanek 2009/8/18: added (String[], double[]) constructor and docs
 */

package edu.uic.ncdm.venn.data;

/** class representing the class combinations and corresponding weights */
public class VennData {
    /** if <code>true</code> then areas is used and <code>data[*][1]</code> is unused, otherwise all areas are implicitly 1 and data contains class pairs */
    public boolean isAreas;
    /** weights for class combinations specified in <code>data[*][0]</code> (if isAreas is <code>true</code>, unused otherwise) */
    public double[] areas;
    /** data is either a list of pairs of class names which all have the weight 1 (isAreas is <code>false</code>) or a list of class combination specifications at [][0] (names are separated by <code>~</code>) with [][1] unused (<code>null</code>, isAreas is <code>true</code>) */
    public String[][] data;

    /** create Venn data by specifying raw contents - no consistency checks are made so it is the responsibility of the caller to make sure the content combination is valid.
     *  @param data (see {@link #data})
     *  @param areas (see {@link #areas})
     *  @param isAreas (see {@limk #isAreas}) */
    public VennData(String[][] data, double[] areas, boolean isAreas) {
        this.data = data;
        this.areas = areas;
        this.isAreas = isAreas;
    }
    
    /** create Venn data from a list of strings specifying the classes and associated areas.
     *  @param data list of string specifying classes/combinations (names must be separated by <code>~</code>)
     *  @param areas areas associated with the classes */
    public VennData(String[] data, double[] areas) {
	this.data = new String[data.length][2];
	for (int i = 0; i < data.length; i++) this.data[i][0] = data[i];
	this.areas = areas;
	this.isAreas = true;
    }

    /** create Venn data from a list of pairs (all will have equal area).
     *  @param pair1 first item of each pair
     *  @param pair2 second item of each pair */
    public VennData(String[] pair1, String[] pair2) {
	this.data = new String[pair1.length][2];
	for (int i = 0; i < data.length; i++) {
	    this.data[i][0] = pair1[i];
	    this.data[i][1] = pair2[i];
	}
	this.isAreas = false;
    }
}
