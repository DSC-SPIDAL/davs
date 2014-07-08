package edu.indiana.soic.spidal.davs;

class Smaller implements Comparable<Smaller>{
    int index;
    double value;

    @Override
    public int compareTo(Smaller o) {
        return Double.compare(value, o.value);
    }
}
