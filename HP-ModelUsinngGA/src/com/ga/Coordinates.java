

package com.ga;

//To help building the plot by storing and creating a data type as Coordinates
public class Coordinates {
    
    private int x;
    private int y;
    
//To give the previous x coordinate
    public int getX() {
        return x;
    }
//To set a new x coordinate
    public void setX(int x) {
        this.x = x;
    }
//To get the previous y coordinate
    public int getY() {
        return y;
    }
//To set a new y coordinate
    public void setY(int y) {
        this.y = y;
    }
    
//To evaluate if the coordinates exist in the path or not   
   @Override
   public boolean equals(Object otherpair){
       
      Coordinates next = (Coordinates) otherpair;
       
       return (this.x == next.x) && (this.y == next.y);
             
   }
}
