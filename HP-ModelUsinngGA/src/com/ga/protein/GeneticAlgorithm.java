package com.ga.protein;
import com.ga.*;
import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

//To draw graphical structures 
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;

public class GeneticAlgorithm extends JComponent {
//initializing the coordinate list

    List<Coordinates> coordinatesList = null; //To maintain a list of coordinates 
    int current = 0;
    static int bestfitness = 0;
    static int finalfitness = 0;
//Initializing Edge with the basic points

    private static class Edge {

        final int x1;
        final int y1;
        final int x2;
        final int y2;
        final Color clr;

        public Edge(int x1, int y1, int x2, int y2, Color clr) {
            this.x1 = x1;
            this.y1 = y1;
            this.x2 = x2;
            this.y2 = y2;
            this.clr = clr;
        }
    }
    //To build a graph with the generated edges
    private final LinkedList<Edge> edges = new LinkedList<>();

    //adding edge with black color
    public void addEdge(int x1, int x2, int x3, int x4) {
        addEdge(x1, x2, x3, x4, Color.black);
    }

    //to add a new edge to the plot with amino acid
    public void addEdge(int x1, int x2, int x3, int x4, Color clr) {
    	 repaint();
        edges.add(new Edge(x1, x2, x3, x4, clr));
        
        repaint();
    }

    //to clear the plot
    public void clearEdges() {
        edges.clear();
        repaint();
    }

    public static void main(String[] args) {
        //Getting the data from Calc class
        Calc.main(args);
        
        // if the coordinate list is ended
        if (!Calc.nextList.isEmpty()) 
        {
        	drawFrame();
        }
    }
        public static void drawFrame()
        {
        	
        	JFrame testFrame = new JFrame();
          //Setting to stop execution once the jframe is closed
          testFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
          //Building the jframe
          final GeneticAlgorithm comp = new GeneticAlgorithm();
          comp.setPreferredSize(new Dimension(1000, 700));
          testFrame.getContentPane().add(comp, BorderLayout.CENTER);
          //setting frame title
          if(Calc.finalIndex == 0)
          {
                 testFrame.setTitle("Protein structure with greatest fitness:- " + (Calc.max));
                 bestfitness = Calc.max;
                 
          }else
          {
        	  testFrame.setTitle("Protein structure with greatest fitness:- " + (Calc.finalIndex));
        	  bestfitness = Calc.finalIndex;
          }
          
       //fitnessArray.add(Calc.finalIndex);
          //initializing the plot
          Calc.previousList = null;
          if (Calc.previousList == null) {
              //getting the first coordinates of the plot
              Calc.previousList = Calc.nextList.get(Calc.mainSubset);
          }
          
          System.out.println("Starting to plot data: " + Calc.previousList.size());
          for (int i = 0; i < (Calc.previousList.size()-1); i++) {
              int x1 = (Calc.previousList.get(i).getX()) * 24;
              int x2 = (Calc.previousList.get(i + 1).getX()) * 24;
              int y1 = (Calc.previousList.get(i).getY()) * 28;
              int y2 = (Calc.previousList.get(i + 1).getY()) * 28;
              System.out.println("(" + x1 + "," + y1 + ")" + "-" + "(" + x2 + "," + y2 + ")");
              Color randomColor = new Color((float) Math.random(), (float) Math.random(), (float) Math.random());
              comp.addEdge(x1, y1, x2, y2, randomColor);
          }
          
          JPanel buttonsPanel = new JPanel();
          testFrame.getContentPane().add(buttonsPanel, BorderLayout.SOUTH);
          testFrame.pack();
          testFrame.setVisible(true);
        }

    
	@Override
    protected void paintComponent(Graphics g1) { //color represention of protein points 
        super.paintComponent(g1);
        Graphics2D g = (Graphics2D) g1.create();
        for (Edge line : edges) {
            g.setColor(Color.BLUE);
            g.setStroke(new BasicStroke(2));
            for (int hpos : Calc.hpPositions) {
                if (current == hpos) {
                    g.setColor(Color.BLACK);
                    break;
                }
            }
            current++;
            g.fillOval(line.x1 + 325, line.y1 + 405, 10, 10);
            g.setColor(Color.BLACK);
            //draw line between two proteins 
            g.drawLine(line.x1 + 325, line.y1 + 405, line.x2 + 325, line.y2 + 405);
            if (current == (Calc.previousList.size() - 1)) {
                for (int hpos : Calc.hpPositions) {
                    if (current == hpos) {
                        g.setColor(Color.ORANGE); //indicating first protein 
                        g.fillOval(line.x2 + 325, line.y2 + 405, 10, 10);
                        break;
                    }
                }
                current = 0;
            }
        }
    }


}

