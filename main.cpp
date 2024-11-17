#include "earclipper.h"

int main (int argc, char** argv) {
    using namespace std;
    if(argc != 2) {
        cerr << "Usage: " << argv[0] << " polygon_csv_filename\n";
        return 1;
    }

    list<Point> points;
    read_from_file(argv[1], points);
    EarClipper clipper(std::move(points));
    clipper();

    return 0;
}
