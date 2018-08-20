package com.ga;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.Scanner;

import com.ga.protein.GeneticAlgorithm;


public class Calc {

    //Initiazing variables for the genetic algorithm
    public static List<Coordinates> previousList = null;
    public List<Integer> previousResult = null;
    public List<Coordinates> previousPoint = null;
    public List<List<Coordinates>> previousPoints = null;
    private int subset1 = 0;
    private int subset2 = 0;
    public static int max=0;
    public static int mainSubset;
    static List<Integer> currentSolutionFitness = null;
    private List<Coordinates> hPoints = null;
    private List<Coordinates> nPoints = null;
    static final int limit = 200000;
    public static int coordinatesCount = 0;
    public static int count = 0;
    public static int nextCount = 1;
    private List<Coordinates> coordinateList = null;
    private List<Coordinates> current = null;
    private static char proteinSequence[] = null;
    private int proteinLength = 0;
    private List<Coordinates> currentCoordinatesList = null;
    private List<Coordinates> currentCoordinate = null;
    static final List<List<Coordinates>> changedList = new ArrayList<>();
    private List<List<Coordinates>> limit1 = null;
    private List<List<Coordinates>> limit2 = null;
    public static int finalIndex = 0;
    public static int finalFitness = 0;
    public static List<List<Coordinates>> nextList = null;
    public static List<Integer> hpPositions = null;
    public static List<List<Coordinates>> changeList = null;
    static private int hCount = 0;
    private int key = 0;
    public static List<Integer> occurances = null;
    static List<Integer> resultList = null;
    private static List<Coordinates> finalCoordinates = null;
    //Get previous coordinates
    public List<Coordinates> getCoordinateList() {
        return coordinateList;
    }
    //Get protein sequence
    public char[] getProtein() {
        return proteinSequence;
    }
    //Setting static values used in plotting the graph  
    static {
        hpPositions = new ArrayList<>();
        occurances = new ArrayList<>();
        changeList = new ArrayList<>();
        finalCoordinates = new ArrayList<>();
        nextList = new ArrayList<>();
        resultList = new ArrayList<>();
        currentSolutionFitness = new ArrayList<>();
    }
    //Construct function
    public Calc() {
        this.proteinLength = proteinSequence.length;
        this.coordinateList = new ArrayList<>();
        this.current = new ArrayList<>();
        this.currentCoordinatesList = new ArrayList<>();
        this.currentCoordinate = new ArrayList<>();
        this.limit1 = new ArrayList<>();
        this.limit2 = new ArrayList<>();
        this.previousResult = new ArrayList<>();
        this.previousPoint = new ArrayList<>();
        this.previousPoints = new ArrayList<>();

    }
    public static void main(String... args) {
        //Starting of the execution
        System.out.println("Enter the protein sequence:");
        Scanner sc = new Scanner(System.in);
        String seq = sc.nextLine();
        //Closing the solution if an empty string is entered
        if (seq.isEmpty()) {
            System.out.println("Empty string cannot be processed");
            System.exit(0);
        } else {
            Calc.proteinSequence = seq.toCharArray();
        }
        
        System.out.println("protein structure length is:  " + proteinSequence.length);
        for (int i = 0; i < proteinSequence.length; i++) {

            if (proteinSequence[i] == 'h') {
                occurances.add(i);
                if (i > 0) {
                    if (proteinSequence[i - 1] == proteinSequence[i]) {
                        hCount++;
                    }
                }
            }
        }

        if (hpPositions.isEmpty()) {
            hpPositions = occurances;
        }
        Calc.generateProteinPoints();

    }
    public static boolean generateProteinPoints(){
    
        int LOWER_BOUND = 10;
        boolean flag = false;
        Calc o = new Calc();
        //Preserving Elite chromozomes
        for (int i = 0; i < LOWER_BOUND; i++) {
            Calc ob = new Calc();
            ob.generateCoordinates();
            ob.mutateStructure();
            if (ob.hPoints != null) {
                o.limit1.add(ob.hPoints);
                o.previousResult.add(o.getStructureFitness(ob.previousPoint));
                o.previousPoints.add(ob.previousPoint);
            }
        }
        //Adding Crossover chromozomes
         max = Collections.max(o.previousResult);
        int indexMax = o.previousResult.indexOf(max);
        o.limit1.remove(indexMax);
        nextList.add(o.previousPoints.get(indexMax));//elite
        o.limit2.add(o.previousPoints.get(indexMax));//coordinates
        Calc.currentSolutionFitness.add(max);//fitness

        int index, init = 0;
        
        for (int step = 0; step < 100; step++) {
            //Generating the best solution
            if (step > 0) {
                List<Coordinates> greats = new ArrayList<>();
                for (Coordinates p : finalCoordinates) {
                    Coordinates temp = new Coordinates();
                    temp.setX(p.getX());
                    temp.setY(p.getY());
                    greats.add(temp);

                }
                Calc.nextList.clear();
                
                Calc.nextList.add(greats);
                Calc.currentSolutionFitness.add(finalIndex);
                o.limit2.add(greats);

            }
            while (flag = (nextList.size() != LOWER_BOUND)) {
                //maximum number of iterations reached
                if (init == limit) {
                    break;
                }
                //Forming new Chromozomes
                
                Calc geno_type = new Calc();
                for (List<Coordinates> lptemp : o.limit1) {
                    List<Coordinates> temp = new ArrayList<>();
                    for (Coordinates p : lptemp) {
                        Coordinates ptemp = new Coordinates();
                        ptemp.setX(p.getX());
                        ptemp.setY(p.getY());
                        temp.add(ptemp);
                    }

                    geno_type.limit1.add(temp);
                }
                //new structure crossovers    
                geno_type.crossoverStructure();
                if (geno_type.nPoints != null) {
                    o.limit2.add(geno_type.nPoints);
                }
                init++;
            }

            nextCount++;
            	GeneticAlgorithm.drawFrame(); //generate graph for each generation 
            	
            
            if (!flag) {
                finalCoordinates.clear();
                Calc ob2 = new Calc();
                //Comparing the fitness of  the new coordinates
                index = ob2.generateBestFit();
                for (Coordinates p : Calc.nextList.get(index)) {
                    Coordinates temp = new Coordinates();
                    temp.setX(p.getX());
                    temp.setY(p.getY());
                    finalCoordinates.add(temp);

                }
                o.limit2.remove(index);

            }

            if (init == limit) {

                init = 0;
                break;
            }
            //Clearing the solution
           // if (step != 0 && step != 199) {
            	if(step!=199) {
                nextList.clear();
                currentSolutionFitness.clear();
                changeList.clear();
                o.previousPoint.clear();
                o.previousResult.clear();
                o.previousPoints.clear();

            }
            o.limit1.clear();
            currentSolutionFitness.clear();
            for (int i = 0; i < 9; i++) {
                Calc ob = new Calc();
                ob.coordinateList = o.limit2.get(i);
                ob.mutateStructure();
                if (ob.hPoints != null) {
                    o.limit1.add(ob.hPoints);
                    o.previousResult.add(o.getStructureFitness(ob.coordinateList));
                    o.previousPoints.add(ob.coordinateList);

                }
            }

            o.limit2.clear();

        }  
      return true;
    }
    //Generate the direction of the next edge   
    public int getStep() {
        Random rotation = new Random();
        
        return rotation.nextInt(4) + 1;
    }
    //Setting the point for the amino acid point    
    public boolean checkPoint(Coordinates p) {
        boolean overlap = false;

        for (Coordinates q : this.coordinateList) {
            if (p.equals(q)) {

                overlap = true;
                break;
            }

        }
        return overlap;
    }
    //Generating the coordinates for the amino acid point
    public boolean generateCoordinates() {

        int previous = 0, current, X1 = 0, Y1 = 0, x = 0, y = 0, i = 0;

        while (i != this.proteinLength) {
            Coordinates p = new Coordinates();
            if (i == 0) {
                p.setX(x);
                p.setY(y);
                this.coordinateList.add(p);
                i++;
            } else {
                current = this.getStep();
                //Checking the SAW condition and obtaining the next direction
                switch (current) {
                    case 1: 
                        x = X1 + 1;
                        y = Y1;
                        p.setX(x);
                        p.setY(y);
                        if (this.checkPoint(p)) {
                            this.coordinateList.clear();
                            previous = 0;
                            X1 = 0;
                            Y1 = 0;
                            x = 0;
                            y = 0;
                            i = 0;
                            break;
                        } else {
                            this.coordinateList.add(p);
                            X1 = x;
                            previous = current;
                            i++;
                            break;
                        }
                    case 2:
                        x = X1 - 1;
                        y = Y1;
                        p.setX(x);
                        p.setY(y);

                        if (this.checkPoint(p)) {
                            this.coordinateList.clear();
                            previous = 0;
                            X1 = 0;
                            Y1 = 0;
                            x = 0;
                            y = 0;
                            i = 0;
                            break;
                        } else {
                            this.coordinateList.add(p);

                            X1 = x;
                            previous = current;
                            i++;
                            break;
                        }

                    case 3: 
                        y = Y1 + 1;
                        x = X1;
                        p.setX(x);
                        p.setY(y);

                        if (this.checkPoint(p)) {
                            this.coordinateList.clear();
                            previous = 0;
                            X1 = 0;
                            Y1 = 0;
                            x = 0;
                            y = 0;
                            i = 0;
                            break;
                        } else {
                            this.coordinateList.add(p);

                            Y1 = y;
                            previous = current;
                            i++;
                            break;
                        }

                    case 4: 
                        y = Y1 - 1;
                        x = X1;
                        p.setX(x);
                        p.setY(y);

                        if (this.checkPoint(p)) {
                            this.coordinateList.clear();
                            previous = 0;
                            X1 = 0;
                            Y1 = 0;
                            x = 0;
                            y = 0;
                            i = 0;
                            break;
                        } else {
                            this.coordinateList.add(p);

                            Y1 = y;
                            previous = current;
                            i++;
                            break;
                        }

                }

            }
        }
        System.out.println("Parent Coordinates# " + coordinatesCount++);
        for (Coordinates p : this.coordinateList) {
            Coordinates temp = new Coordinates();
            temp.setX(p.getX());
            temp.setY(p.getY());
            this.previousPoint.add(temp);
        }
        return true;

    }
    //Calculating the plot fitness
    public int getStructureFitness(List<Coordinates> list) {

        int fit = 0;
        if (occurances.isEmpty() || list == null) {
            return 0;
        } else {
            //listing each hydrophobic and hydrophilic occurance 
            for (int hp : occurances) {
                Coordinates originpair = list.get(hp);
                //Listing the hidrophobic occurance
                for (int h : occurances) {
                    Coordinates nextpair = list.get(h);
                    int x1 = originpair.getX();
                    int y1 = originpair.getY();

                    int x2 = nextpair.getX();
                    int y2 = nextpair.getY();

                    double dsq = (Math.pow((x1 - x2), 2) + Math.pow((y1 - y2), 2));
                    float distance = (float) Math.sqrt(dsq);

                    if (distance == 1.0) {
                        fit++;
                    }

                }
            }
            //Calculating Fitness
            fit = (fit / 2) - hCount;
            resultList.add(fit);
            finalFitness = fit;

            return fit;
        }
    }
    //Get fitness after crossover
    public int getCrossoverResult(List<Coordinates> list) {

        int fit = 0;
        if (occurances.isEmpty() || list == null) {
            return 0;
        } else {
            for (int hp : occurances) {
                Coordinates originpair = list.get(hp);
                for (int h : occurances) {
                    Coordinates nextpair = list.get(h);
                    int x1 = originpair.getX();
                    int y1 = originpair.getY();

                    int x2 = nextpair.getX();
                    int y2 = nextpair.getY();
                    //Calculating distace between the pairs
                    double dsq = (Math.pow((x1 - x2), 2) + Math.pow((y1 - y2), 2));
                    float distance = (float) Math.sqrt(dsq);
                    if (distance == 1.0) {
                        fit++;
                    }
                    
                }
            }
            //Calculating fitness
            fit = (fit / 2) - hCount;
            finalFitness = fit;
            return fit;
        }
    }
    //To get fitness after mutation of the structure
    public boolean mutateStructure() {
        System.out.println("Mutation Started for Generation: " + nextCount);
        this.generateRotation();
        return true;
    }
    //Get the structure after crossover 
    public boolean crossoverStructure() {
        count = 1;
        int pivot1 = this.generatePivotForCrossover();
        int pivot2 = pivot1 + 1;

        subset1 = this.generateIndex();
        subset2 = this.generateIndex();
        //Cross joining two solutions at a pivot point
        List<Coordinates> list1 = this.limit1.get(subset1).subList(0, pivot1 + 1);
        List<Coordinates> list2 = this.limit1.get(subset2).subList(pivot2, this.proteinLength);
        
        Coordinates o = new Coordinates();
        o.setX(this.limit1.get(subset1).get(pivot2).getX());
        o.setY(this.limit1.get(subset1).get(pivot2).getY());
        //Get the previous coordinate points
        int offset_x = list2.get(0).getX() - o.getX();
        int offset_y = list2.get(0).getY() - o.getY();
        //Saving newly formed chromozome coordinates
        for (Coordinates z : list2) {
            int x = z.getX() - offset_x;
            int y = z.getY() - offset_y;
            z.setX(x);
            z.setY(y);

        }
        //Joining two parts of different chromozomes
        list1.addAll(list2);
        
        this.currentCoordinatesList.addAll(list1);

        this.key = pivot2;
        //Checking if the crossover meets the restricts of genetic algrithm
        if (this.rotateCrossoverSAW(0)) {
        	//GeneticAlgorithm.drawFrame();
        	System.out.println("Crossover Success");
            return true;

        } else {
            return false;
        }

    }
    //To check the new crossover chromozome if it meets the restrictions of genetic algorithm
    public boolean rotateCrossoverSAW(int step) {

        count++;
        Coordinates origin = this.currentCoordinatesList.get(this.key);

        int X, Y;
        boolean flag = false;
        for (Coordinates p : this.currentCoordinatesList) {
            Coordinates temp = new Coordinates();
            temp.setX(p.getX());
            temp.setY(p.getY());
            this.currentCoordinate.add(temp);
        }

        if (currentCoordinate.isEmpty()) {
            System.out.println("Error:Holder is empty!");

        } else if (step == 0) {

            if (flag = !this.validateCoordinatestoAdd(currentCoordinate)) {
                currentCoordinate.clear();
                return this.rotateCrossoverSAW(1);

            }

        } else if (step == 1) {
            for (int i = this.key + 1; i < proteinSequence.length; i++) {
                Coordinates p = new Coordinates();
                X = (this.currentCoordinatesList.get(i).getY() * -1) + (origin.getY() + origin.getX());
                Y = (this.currentCoordinatesList.get(i).getX()) + (origin.getY() - origin.getX());
                p.setX(X);
                p.setY(Y);
                this.currentCoordinate.set(i, p);
            }
            if (flag = !this.validateCoordinatestoAdd(currentCoordinate)) {
                currentCoordinate.clear();
                return this.rotateCrossoverSAW(3);

            }

        } else if (step == 3) {

            for (int i = this.key + 1; i < proteinSequence.length; i++) {
                Coordinates p = new Coordinates();
                X = (this.currentCoordinatesList.get(i).getY()) + (origin.getX() - origin.getY());
                Y = (this.currentCoordinatesList.get(i).getX() * -1) + (origin.getX() + origin.getY());
                p.setX(X);
                p.setY(Y);
                this.currentCoordinate.set(i, p);
            }
            if (flag = !this.validateCoordinatestoAdd(currentCoordinate)) {
                currentCoordinate.clear();
            }
        }
        return !flag;
    }
    //To checck the coordinate can be added to the plot
    public boolean validateCoordinatestoAdd(List<Coordinates> current) {
        int X = 0, Y = 0, fp1, fp2, newgen;
        double avg;
        boolean flag = false;
        
        if (this.VerifyCrossoverCoordinatesOverlap() || this.VerifyCrossoverCoordinatesDistance()) {
            current.clear();
            return flag;
        } else {
            fp1 = this.getStructureFitness(changeList.get(subset1));
            fp2 = this.getStructureFitness(changeList.get(subset2));
            avg = (fp1 + fp2) / 2;
            if ((newgen = this.getCrossoverResult(current)) > avg) {
                nextList.add(current);
                nPoints = current;
                currentSolutionFitness.add(newgen);
                flag = true;
            } else {
                current.clear();
                flag = false;
            }
            return flag;
        }
    }
    //To Rotate the plot
    public boolean generateRotation() {
        for (Coordinates p : this.coordinateList) {
            Coordinates temp = new Coordinates();
            temp.setX(p.getX());
            temp.setY(p.getY());
            this.current.add(temp);
        }
        
      this.key = this.generatePivot();
        
        Coordinates origin = this.coordinateList.get(this.key);
        int X, Y, f1, f2;
        boolean flag = false;
        //To rotate the plot in given angles
        switch (this.generateAngle() * 90) {
            case 90:
                for (int i = this.key + 1; i < proteinSequence.length; i++) {
                    Coordinates p = new Coordinates();
                    X = (this.coordinateList.get(i).getY() * -1) + (origin.getY() + origin.getX());
                    Y = (this.coordinateList.get(i).getX()) + (origin.getY() - origin.getX());
                    p.setX(X);
                    p.setY(Y);
                    this.current.set(i, p);
                }
                if (current.isEmpty()) {
                    System.out.println("Error:Holder is empty!");
                    break;
                } else {
                    if (this.checkOverlap() || this.VerifyCoordinatesDistance()) {
                        current.clear();
                        return this.generateRotation();
                    } else {
                        flag = true;
                        if ((f1 = this.getStructureFitness(this.coordinateList)) <= (f2 = this.getStructureFitness(this.current))) {
                            hPoints = this.current;
                            changeList.add(this.current);
                            resultList.add(f2);
                        } else {
                            hPoints = this.coordinateList;
                            changeList.add(this.coordinateList);
                            resultList.add(f1);
                        }
                    }
                    break;
                }
            case 180:
                for (int i = this.key + 1; i < proteinSequence.length; i++) {
                    Coordinates p = new Coordinates();
                    X = (this.coordinateList.get(i).getX() * -1) + (2 * origin.getX());
                    Y = (this.coordinateList.get(i).getY() * -1) + (2 * origin.getY());
                    p.setX(X);
                    p.setY(Y);
                    this.current.set(i, p);
                }
                if (current.isEmpty()) {
                    System.out.println("Error:Holder is empty!");
                    break;
                } else {
                    if (this.checkOverlap() || this.VerifyCoordinatesDistance()) {
                        current.clear();
                        return this.generateRotation();
                    } else {
                        flag = true;
                        if ((f1 = this.getStructureFitness(this.coordinateList)) <= (f2 = this.getStructureFitness(this.current))) {
                            hPoints = this.current;
                            changeList.add(this.current);
                            resultList.add(f2);
                        } else {
                            hPoints = this.coordinateList;
                            changeList.add(this.coordinateList);
                            resultList.add(f1);
                        }
                    }
                    break;
                }
            case 270:
                for (int i = this.key + 1; i < proteinSequence.length; i++) {
                    Coordinates p = new Coordinates();
                    Coordinates q = new Coordinates();
                    X = (this.coordinateList.get(i).getY()) + (origin.getX() - origin.getY());
                    Y = (this.coordinateList.get(i).getX() * -1) + (origin.getX() + origin.getY());
                    p.setX(X);
                    p.setY(Y);
                    this.current.set(i, p);
                }
                if (current.isEmpty()) {
                    System.out.println("Error:Holder is empty!");
                    break;
                } else {
                    if (this.checkOverlap() || this.VerifyCoordinatesDistance()) {
                        current.clear();
                        return this.generateRotation();
                    } else {
                        flag = true;
                        if ((f1 = this.getStructureFitness(this.coordinateList)) <= (f2 = this.getStructureFitness(this.current))) {
                            hPoints = this.current;
                            changeList.add(this.current);
                            resultList.add(f2);
                        } else {
                            hPoints = this.coordinateList;
                            changeList.add(this.coordinateList);
                            resultList.add(f1);
                        }
                    }
                    break;
                }
            default:
                System.out.println("Angle Does Not Exist");

        }

        return flag;
    }
    //Generating anglel for the chromozome
    public int generateAngle() {

        Random r = new Random();
        int angle = r.nextInt(3) + 1;
        return angle;
    }
    
    public int generatePivot() {

        Random r = new Random();

        return r.nextInt(proteinLength - 1);

    }
    public int generateAngleForCrossover() {

        Random r = new Random();

        return ((2 * r.nextInt(3 - 1)) + 1);
    }
    //Generating Pivot point for cross over funtion
    public int generatePivotForCrossover() {
        Random r = new Random();

        return r.nextInt(proteinLength - 1);
    }
   
    //To Veify the overlap while building the Plot
    public boolean checkOverlap() {
        boolean flag = false;
        for (int i = key + 1; i < proteinSequence.length; i++) {
            for (int step = 0; step < key + 1; step++) {
                if (current.get(i).equals(current.get(step))) {
                    flag = true;
                    break;
                }
            }
            if (flag) {
                break;
            }
        }

        return flag;

    }
    //To verify if the coordidnates overlap after the crossover
    public boolean VerifyCrossoverCoordinatesOverlap() {
        boolean flag = false;
        for (int i = key + 1; i < proteinSequence.length; i++) {
            for (int step = 0; step < key + 1; step++) {
                if (currentCoordinate.get(i).equals(currentCoordinate.get(step))) {
                    flag = true;
                    break;
                }
            }
            if (flag) {
                break;
            }
        }

        return flag;

    }
    //To verify the coordinte mismatch
    public boolean VerifyCoordinatesDistance() {

        Coordinates p, q;
        int remaining_length = (proteinLength - key - 2);
        double dst;
        float distance, dsum = 0;

        for (int i = key + 1; i < proteinSequence.length - 1; i++) {
            p = current.get(i);
            q = current.get(i + 1);
            dst = (Math.pow((p.getX() - q.getX()), 2) + Math.pow((p.getY() - q.getY()), 2));
            distance = (float) Math.sqrt(dst);

            dsum = dsum + distance;
        }
        return !(dsum == (float) remaining_length);

    }
    //To check the mismatched coordinates after cross over
    

    
    //To create random index
    public int generateIndex() {

        Random r = new Random();
        if (count == 1) {
            return r.nextInt(9);
        } else {
            return r.nextInt(10);
        }
    }
    //To Get the best solution
    public int generateBestFit() {
        finalIndex = Collections.max(currentSolutionFitness);
        mainSubset = currentSolutionFitness.indexOf(finalIndex);
        return mainSubset;

    }

public boolean VerifyCrossoverCoordinatesDistance() {

        Coordinates p, q;
        int remaining_length = (proteinLength - key);

        double dst;
        float distance, dsum = 0;

        for (int i = key - 1; i < proteinSequence.length - 1; i++) {
            p = currentCoordinate.get(i);
            q = currentCoordinate.get(i + 1);
            dst = (Math.pow((p.getX() - q.getX()), 2) + Math.pow((p.getY() - q.getY()), 2));
            distance = (float) Math.sqrt(dst);
            dsum = dsum + distance;
        }
        return !(dsum == (float) remaining_length);

    }
}