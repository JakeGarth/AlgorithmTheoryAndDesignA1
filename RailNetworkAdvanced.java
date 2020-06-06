package assg1p2;

import assg1p1.Station;

import java.io.*;
import java.util.*;


public class RailNetworkAdvanced {

    private TreeMap<String, Station> stationList;
    public HashSet<String> unsettled;
    public HashMap<Station, Integer> distance = new HashMap<Station, Integer>();
    public HashMap<String, String> tracker = new HashMap<String, String>();

    public HashMap<String, String[]> lineData = new HashMap<>();
    public HashMap<String, String[]> lineStationData = new HashMap<>();
    public String[] firstRow = new String[20];


    public RailNetworkAdvanced(String trainData, String connectionData, String lineData) {
        stationList = new TreeMap<>();


        try {
            readLinesData(lineData);
            readStationData(trainData);
            readConnectionData(connectionData);
        } catch (IOException e) {
            System.out.println("Exception encountered: " + e);
        }
    }

    public void reset() {
        unsettled = new HashSet<String>();
        tracker = new HashMap<String, String>();
        distance = new HashMap<Station, Integer>();

    }

    /**
     * Reads the CSV file containing information about the lines
     *
     * @param infile
     * @throws IOException
     */
    public void readLinesData(String infile) throws IOException {
        BufferedReader in = new BufferedReader(new FileReader(infile));

        String s;

        s = in.readLine();

        while ((s = in.readLine()) != null) {

            String[] data = s.split(",", -1);
            //puts in all the data from the "lines_data.csv" with the "code" as the key
            lineData.put(data[0], data);

        }
        in.close();
    }


    /**
     * Reads the CSV file containing information about the stations and
     * populate the TreeMap<String,Station> stationList. Each row of
     * the CSV file contains the name, latitude and longitude coordinates
     * of the station.
     * <p>
     * You need to make a Station object for each row and add it to the
     * TreeMap<String,Station> stationList where the key value is the
     * name of the station (as a String).
     *
     * @param infile the path to the file
     * @throws IOException if the file is not found
     */
    public void readStationData(String infile) throws IOException {
        BufferedReader in = new BufferedReader(new FileReader(infile));
        /*
         * INSERT YOUR CODE HERE
         */
        String s;
        s = in.readLine();

        String[] firstLine = s.split(",", -1);
        firstRow = firstLine;

        while ((s = in.readLine()) != null) {

            String[] data = s.split(",", -1);

            lineStationData.put(data[0], data);
            Station stat = new Station(data[0], Double.valueOf(data[1]), Double.valueOf(data[2]));
            stationList.put(data[0], stat);

        }
        in.close();
    }

    /**
     * Reads the CSV file containing information about connectivity between
     * adjacent stations, and update the stations in stationList so that each
     * Station object has a list of adjacent stations.
     * <p>
     * Each row contains two Strings separated by comma. To obtain the distance
     * between the two stations, you need to use the latitude and longitude
     * coordinates together with the computeDistance() methods defined below
     *
     * @param infile the path to the file
     * @throws IOException if the file is not found
     */
    public void readConnectionData(String infile) throws IOException {
        BufferedReader in = new BufferedReader(new FileReader(infile));
        /*
         * INSERT YOUR CODE HERE
         */
        String s;
        while ((s = in.readLine()) != null) {
            String[] data = s.split(",", -1);

            Station station1 = stationList.get(data[0]);
            Station station2 = stationList.get(data[1]);

            int distance = computeDistance(station1.getLatitude(), station1.getLongitude(),
                    station2.getLatitude(), station2.getLongitude());

            station1.addNeighbour(station2, distance);
            station2.addNeighbour(station1, distance);

            // put only after we have taken old station and updated it
            stationList.put(data[0], station1);
            stationList.put(data[1], station2);
        }
        in.close();
    }

    /**
     * Given the latitude and longitude coordinates of two locations x and y,
     * return the distance between x and y in metres using Haversine formula,
     * rounded down to the nearest integer.
     * <p>
     * Note that two more methods are provided below for your convenience
     * and you should not directly call this method
     * <p>
     * source://www.geeksforgeeks.org/haversine-formula-to-find-distance-between-two-points-on-a-sphere/
     *
     * @param lat1 latitude coordinate of x
     * @param lon1 longitude coordinate of x
     * @param lat2 latitude coordinate of y
     * @param lon2 longitude coordinate of y
     * @return distance betwee
     */
    public static int computeDistance(double lat1, double lon1, double lat2, double lon2) {
        // distance between latitudes and longitudes
        double dLat = Math.toRadians(lat2 - lat1);
        double dLon = Math.toRadians(lon2 - lon1);

        // convert to radians
        lat1 = Math.toRadians(lat1);
        lat2 = Math.toRadians(lat2);

        // apply formulae
        double a = Math.pow(Math.sin(dLat / 2), 2) +
                Math.pow(Math.sin(dLon / 2), 2) *
                        Math.cos(lat1) *
                        Math.cos(lat2);
        double rad = 6371.0;
        Double c = 2 * Math.asin(Math.sqrt(a));
        Double distance = rad * c * 1000;
        return distance.intValue();
    }

    /**
     * Compute the distance between two stations in metres, where the stations
     * are given as String objects
     *
     * @param a the first station
     * @param b the second station
     * @return the distance between the two stations in metres
     */
    public int computeDistance(String a, String b) {
        Station u = stationList.get(a);
        Station v = stationList.get(b);
        return computeDistance(u.getLatitude(), u.getLongitude(),
                v.getLatitude(), v.getLongitude());
    }

    /**
     * Compute the distance between two stations in metres, where the stations
     * are given as Station objects
     *
     * @param a the first station
     * @param b the second station
     * @return the distance between the two stations in metres
     */
    public int computeDistance(Station a, Station b) {
        return computeDistance(a.getLatitude(), a.getLongitude(),
                b.getLatitude(), b.getLongitude());
    }

    /**
     * The method finds the shortest route (in terms of distance travelled)
     * between the origin station and the destination station.
     * The route is returned as an ArrayList<String> containing the names of
     * the stations along the route, including the origin and the destination
     * stations.
     * <p>
     * If the route cannot be completed (there is no path between origin and
     * destination), then return an empty ArrayList<String>
     * <p>
     * If the origin or the destination stations are not in the list of stations,
     * return an empty ArrayList<String>.
     * <p>
     * If the origin and the destination stations are the same, return an
     * ArrayList<String> containing the station.
     *
     * @param origin      the starting station
     * @param destination the destination station
     * @return
     */
    //THIS IS ROUTE MIN DISTANCE
    public ArrayList<String> routeMinDistance(String origin, String destination) {
        reset();
        unsettled = setAllUnmarked();
        setAllDistance(stationList.get(origin));
        makingDistances(origin, destination);
        if (!stationList.containsKey(origin) || !stationList.containsKey(destination)) {
            return new ArrayList<String>();
        }
        if (origin.equals(destination)) {
            ArrayList<String> ans = new ArrayList<String>();
            ans.add(origin);
            return ans;
        }

        ArrayList<String> answer = trackerMaker(tracker, destination, origin);
        Collections.reverse(answer);
        return answer;
    }


    public void setAllDistance(Station origin) {
        for (Map.Entry<String, Station> entry : stationList.entrySet()) {
            if (entry.getValue() == origin) {
                distance.put(entry.getValue(), 0);
            } else {
                distance.put(entry.getValue(), Integer.MAX_VALUE);
            }

        }
    }

    public HashSet<String> setAllUnmarked() {
        HashSet<String> test = new HashSet<String>();
        for (Map.Entry<String, Station> entry : stationList.entrySet()) {
            test.add(entry.getKey());
        }
        return test;
    }

    public void makingDistances(String origin, String destination) {

        unsettled.remove(origin);

        if (unsettled.contains(destination)) {
            TreeMap<Station, Integer> adjStations = stationList.get(origin).getAdjacentStations();
            for (Map.Entry<Station, Integer> entry : adjStations.entrySet()) {
                if (unsettled.contains(entry.getKey().getName())) {

                    if (!tracker.containsKey(entry.getKey().getName())) {
                        tracker.put(entry.getKey().getName(), origin);
                    }
                }
                if (distance.get(entry.getKey()) > (entry.getValue() + distance.get(stationList.get(origin)))) {
                    distance.put(entry.getKey(), entry.getValue() + distance.get(stationList.get(origin)));
                }
            }
            unsettled.remove(origin);
            int smallest = Integer.MAX_VALUE;
            String nameSmallest = "";
            for (String dist : unsettled) {
                if (distance.get(stationList.get(dist)) < smallest) {
                    smallest = distance.get(stationList.get(dist));
                    nameSmallest = stationList.get(dist).getName();
                }
            }

            if (!nameSmallest.equals("")) {
                makingDistances(nameSmallest, destination);
            }
        }
    }

    public ArrayList<String> trackerMaker(HashMap<String, String> map, String destination, String origin) {
        ArrayList<String> test = new ArrayList<>();
        if (destination.equals(origin)) {

            test.add(origin);
            return test;
        }
        test.add(destination);
        if (!map.containsKey(destination)) {
            return new ArrayList<>();
        }
        test.addAll(trackerMaker(map, map.get(destination), origin));
        return test;


    }


    /**
     * Given a route between two stations, compute the total distance
     * of this route.
     *
     * @param path the list of stations in the route (as String objects)
     * @return the length of the route between the first station
     * and last station in the list
     */
    public int findTotalDistance(ArrayList<String> path) {
        int dista = 0;
        if (path == null || path.isEmpty()) {
            return 0;
        }
        int size = path.size();
        for (int i = 0; i < size - 1; i++) {
            int section = computeDistance(path.get(i), path.get(i + 1));
            dista = dista + section;
//			System.out.println("section " + i + ": " + section);
        }
        return dista;
    }


    /**
     * Return the ratio between the length of the shortest route between two
     * stations (in terms of distance) and the actual distance between the
     * two stations (computed using computeDistance())
     * <p>
     * In other words,
     * let d1 = distance of shortest route between the two stations as computed
     * by the method routeMinStop() (from Stage 1).
     * let d2 = distance between two stations as computed by the method
     * computeDistance()
     * <p>
     * The method returns d1/d2 (as a double)
     *
     * @param origin      the starting station
     * @param destination the ending station
     * @return s            the ratio d1/d2 as explained above
     */
    public double computeRatio(String origin, String destination) {
        Station orig = stationList.get(origin);
        Station dest = stationList.get(destination);
        int currentd = findTotalDistance(routeMinDistance(origin, destination));
        int newd = computeDistance(orig.getLatitude(), orig.getLongitude(), dest.getLatitude(), dest.getLongitude());
        return (float) currentd / newd;
    }


    /**
     * Return the ratio as computed by computeRatio() method for all
     * pairs of station in the rail network.
     * <p>
     * The ratios should be stored in a HashMap<String,HashMap<String,Double>>,
     * that is, the ratio between station a and b can be obtained by calling
     * <p>
     * computeAllRatio().get(a).get(b)
     *
     * @return a hashmap containing the ratios
     */
    public HashMap<String, HashMap<String, Double>> computeAllRatio() {

        // the object we are returning
        HashMap<String, HashMap<String, Double>> a = new HashMap<>();


        // put all the 'a keys' in the hashmap before we put the data in, i think it makes it easier this way
        // also find stations with smalllest lat and long
        String lowlat = stationList.firstEntry().getValue().getName();
        String highlat = stationList.firstEntry().getValue().getName();
        String lowlong = stationList.firstEntry().getValue().getName();
        String highlong = stationList.firstEntry().getValue().getName();
        ArrayList<String> mystations1 = new ArrayList<>();
        ArrayList<String> mystations2 = new ArrayList<>();

        for (Map.Entry<String, Station> station : stationList.entrySet()) {
            String stationName = station.getKey();
            a.put(stationName, new HashMap<>());

            if (stationList.get(lowlat).getLatitude() > station.getValue().getLatitude()) {
                lowlat = stationName;
            }
            if (stationList.get(highlat).getLatitude() < station.getValue().getLatitude()) {
                highlat = stationName;
            }
            if (stationList.get(lowlong).getLongitude() > station.getValue().getLongitude()) {
                lowlong = stationName;
            }
            if (stationList.get(highlong).getLongitude() < station.getValue().getLongitude()) {
                highlong = stationName;
            }
            mystations2.add(stationName);
            mystations1.add(stationName);
        }

        // put the opposing corner stations at the beginning of both the ArrayLists
        // the goal of this is to reduce computational redundancy by doing the longer routes first
        mystations1.add(0, lowlat);
        mystations1.add(0, lowlong);
        mystations2.add(0, highlat);
        mystations2.add(0, highlong);

        for (int i = 0; i < mystations1.size(); i++) {
            for (int j = 0; j < mystations2.size(); j++) {
                String start = mystations1.get(i);
                String stop = mystations2.get(j);
                // don't do any of this work if it has already been done
                if (!a.get(start).containsKey(stop)) {
                    ArrayList<String> prefill = routeMinDistance(start, stop);
                    Integer size = prefill.size();
                    // 2D array of distances
                    int[][] distances = new int[size][size];
                    for (int col = 0; col < size; col++) {
                        for (int row = 0; row < size; row++) {
                            if (row > col) {
                                //check if the pair has already been done otherwise skip
                                //this is the first one of the column
                                if (row == col + 1) {
                                    distances[row][col] = computeDistance(prefill.get(row), prefill.get(col));
                                } else {
                                    distances[row][col] = distances[row - 1][col] + computeDistance(prefill.get(row), prefill.get(row - 1));
                                }
                                //calculate ratio and put it in the hashmap
                                // lets make this easier. everything gets a veriable
                                String st1 = prefill.get(col);
                                String st2 = prefill.get(row);
                                Integer railDist = distances[row][col];
                                Integer airDist = computeDistance(st1, st2);
                                //calculate ratio
                                Double ratio = railDist / (double) airDist;
                                // put it in the hashmap
                                a.get(st1).put(st2, ratio);
                                a.get(st2).put(st1, ratio);
                            }
                        }
                    }

                }
            }
        }

        return a;
    }


    //THIS IS MIN STOPS

    public void setAllDistanceStop(Station origin) {
        for (Map.Entry<String, Station> entry : stationList.entrySet()) {
            if (entry.getValue() == origin) {
                distance.put(entry.getValue(), 0);
            } else {
                distance.put(entry.getValue(), Integer.MAX_VALUE);
            }

        }
    }

    public ArrayList<String> trackerMakerStop(HashMap<String, String> map, String destination, String origin) {
        ArrayList<String> test = new ArrayList<>();
        if (destination.equals(origin)) {
            test.add(origin);
            return test;
        }
        if (!map.containsKey(destination)) {
            return new ArrayList<>();
        }
        test.add(destination);

        test.addAll(trackerMakerStop(map, map.get(destination), origin));

        return test;

    }

    public ArrayList<String> routeMinStop(String origin, String destination) {
        unsettled = setAllUnmarked();
        setAllDistanceStop(stationList.get(origin));
        makingDistancesStop(origin, destination);
        if (!stationList.containsKey(origin) || !stationList.containsKey(destination)) {
            return new ArrayList<String>();
        }
        if (origin.equals(destination)) {
            ArrayList<String> ans = new ArrayList<String>();
            ans.add(origin);
            return ans;
        }

        ArrayList<String> answer = trackerMakerStop(tracker, destination, origin);
        Collections.reverse(answer);
        return answer;
    }

    public void makingDistancesStop(String origin, String destination) {

        unsettled.remove(origin);
        if (unsettled.contains(destination)) {
            TreeMap<Station, Integer> adjStations = stationList.get(origin).getAdjacentStations();

            for (Map.Entry<Station, Integer> entry : adjStations.entrySet()) {
                if (unsettled.contains(entry.getKey().getName())) {
                    if (!tracker.containsKey(entry.getKey().getName())) {
                        tracker.put(entry.getKey().getName(), origin);
                    }
                }

                if (distance.get(entry.getKey()) > (1 + distance.get(stationList.get(origin)))) { //exact same thing, just changed it to +1 for each stop
                    distance.put(entry.getKey(), 1 + distance.get(stationList.get(origin)));
                }
            }
            unsettled.remove(origin);
            int smallest = Integer.MAX_VALUE;
            String nameSmallest = "";
            for (String dist : unsettled) {
                if (distance.get(stationList.get(dist)) < smallest) {
                    smallest = distance.get(stationList.get(dist));
                    nameSmallest = stationList.get(dist).getName();
                }
            }

            if (!nameSmallest.equals("")) {
                makingDistancesStop(nameSmallest, destination);
            }
        }
    }


    /**
     * The method finds the shortest route (in terms of number of stops)
     * between the origin station and the destination station, taking
     * into account the available routes in the rail network.
     * <p>
     * The route is returned as an ArrayList<String> containing the lines taken,
     * any transfer between lines, and the names of the stations on each line,
     * including the origin and the destination stations.
     * <p>
     * Please see the assignment specification for more information.
     * <p>
     * If the route cannot be completed (there is no path between origin and
     * destination), then return an empty ArrayList<String>
     * <p>
     * If the origin or the destination stations are not in the list of stations,
     * return an empty ArrayList<String>.
     * <p>
     * If the origin and the destination stations are the same, return an
     * ArrayList<String> containing the station.
     *
     * @param origin      the starting station
     * @param destination the end station
     * @return the route taken
     */
    public ArrayList<String> routeMinStopWithRoutes(String origin, String destination) {
        if(!stationList.containsKey(origin)||!stationList.containsKey(destination)){
            return new ArrayList<>();
        }
        reset(); //clears all the global data
        ArrayList<String> stops = routeMinStop(origin, destination); //Uses Jake Garth's routeMinStop from Stage 1,
        //This gives an ArrayList of the min stop route that will be taken
        ArrayList<String> first = new ArrayList<>();
        first.addAll(stops);
        first = LineList(first); //Apply the "LineList()" method to the original ArrayList. This will include the lines taken.
        return first;
    }

    public String findLine(String stationName, String nextStation) {
        int size = firstRow.length - 1;
        String line = null;
        for (int i = 3; i < size; i++) {
            //lineStationData stores the lines on which a Station exists on
            //if two adjacent Stations run on the same line, this line can be used
            if (!lineStationData.get(stationName)[i].equals("") && !lineStationData.get(nextStation)[i].equals("")) {
                //"line" stores the line being used, by using the lineData global data structure
                //this "if" statement determines the direction in which the route is taking, based on how many stops away
                //the Stations are from the origin of the train Line
                if (Integer.parseInt(lineStationData.get(stationName)[i]) < (Integer.parseInt(lineStationData.get(nextStation)[i]))) {
                    line = lineData.get(firstRow[i])[1] + " towards " + lineData.get(firstRow[i])[3] + " from " + lineData.get(firstRow[i])[2];
                } else {
                    line = lineData.get(firstRow[i])[1] + " towards " + lineData.get(firstRow[i])[2] + " from " + lineData.get(firstRow[i])[3];
                }
            }
        }
        return line;
    }

    public ArrayList<String> LineList(ArrayList<String> original) {
        String lastLineRecorded = findRecentLine(original);  //findRecentLine finds the most recent line recorded in the ArrayList
        int i = 0;
        while (i < original.size() - 1) {  //iterates through the entire ArrayList
            if (original.get(i).contains("towards")) {  //if the current String is a line transfer, skip this
                i++; //findLine finds the line on which two adjacent stations share
                //if the line of the current two stations is different to the last recorded line, insert the new line in the ArrayList
            } else if (!findLine(original.get(i), original.get(i + 1)).equals(lastLineRecorded)) {
                List<String> first = original.subList(0, i);
                if (i > 0) {
                    first = original.subList(0, i + 1);
                }
                //"first" is the first part of the ArrayList
                //"second" is the second part of the ArrayList
                List<String> second = original.subList(i, original.size());
                ArrayList<String> firstArray = new ArrayList<>(first);
                ArrayList<String> secondArray = new ArrayList<>(second);
                //Converting the Lists to ArrayLists
                String newLine = findLine(original.get(i), original.get(i + 1));
                //newLine is the line and direction of travel
                original.clear();
                original.addAll(firstArray);
                original.add(newLine);
                original.addAll(secondArray);
                //the ArrayList has now been correctly merged together
                lastLineRecorded = newLine;
                //Updated the "lastLineRecorded" to check for changes in the track
                i += 2;
                //Skip two positions as the length of the ArrayList just grew by one

            } else {
                i++;
                //if no update to the train line, just iterate to next station

            }
            System.out.println(original);
        }
        //when we reach the end of the ArrayList, return it.
        return original;
    }

    public String findRecentLine(ArrayList<String> original) {
        //simply iterate through the ArrayList and finds the last Line recorded
        //checks if "towards" is in the String to determine if its a Station or a Line
        String s = "";
        for (int i = 0; i < original.size(); i++) {
            if (original.get(i).contains("towards")) {
                s = original.get(i);
            }
        }
        return s;
    }

}