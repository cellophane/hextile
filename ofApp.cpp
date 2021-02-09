#include "ofApp.h"
struct hexInfo {
	int squaredDistance;
	int pieceID;
};
struct CubeCoord {
	int x{0};
	int y{0};
	int z{0};

};
bool operator<(const CubeCoord& l, const CubeCoord& r) {
	return (l.x < r.x || (l.x == r.x && l.y < r.y));
}
struct Point {
	float x;
	float y;
};
const float pi = 3.14159;
//radius of the major hexagon
int hexSize = 100;

//minor hex size in pixels (center to vertex)
int hexSizePix = 1;
float hexWidth = sqrt(3) * hexSizePix;
float hexHeight = 2 * hexSizePix;

//cube to cartesian coordinate vectors
float xCube[2] = { hexWidth / 2.,3. / 4. * hexHeight };
float yCube[2] = { -hexWidth / 2.,3. / 4. * hexHeight };

//pixel coordinate of center of major hexagon
ofPoint pixelCenter((hexSize+.5)* hexWidth, hexSize* hexHeight);
CubeCoord center = { 0,0,0 };
//cardinal hex directions
const CubeCoord cubeDirections[6] = { {1,-1,0},{1,0,-1},{0,1,-1},{-1,1,0},{-1,0,1},{0,-1,1} };

//holds the centers of the adjacent major hexagons
vector<CubeCoord> mirrorCenters;
//hexagons to pieces and squared distances
map<CubeCoord, pair<int, int>> hexes;
//piece origins
vector<CubeCoord> pieceOrigins;

vector<ofColor> colors;

CubeCoord rotateClockwise(CubeCoord a) {
	return CubeCoord{ -a.z,-a.x,-a.y };
}
CubeCoord cubeAdd(CubeCoord a, CubeCoord b) {
	return { a.x + b.x, a.y + b.y, a.z + b.z };
}
CubeCoord cubeSub(CubeCoord a, CubeCoord b) {
	return { a.x - b.x, a.y - b.y, a.z - b.z };
}
int cubeDistance(CubeCoord a, CubeCoord b) {
	return max(max(abs(a.x - b.x), abs(a.y - b.y)), abs(a.z - b.z));
}
int squareCubeEuclideanDistance(CubeCoord a, CubeCoord b) {
	int d = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) - (a.x - b.x) * (a.y - b.y);
	if (d >= hexSize * hexSize) {
		for (auto c : mirrorCenters) {
			CubeCoord a0 = cubeAdd(a, c);
			d = min(d, (a0.x - b.x) * (a0.x - b.x) + (a0.y - b.y) * (a0.y - b.y) - (a0.x - b.x) * (a0.y - b.y));
		}
	}

	return d;
}
Point euclideanCenter(CubeCoord a) {
	Point p;
	p.x = a.x * xCube[0] + a.y * yCube[0];
	p.y = a.x * xCube[1] + a.y * yCube[1];
	return p;
}
float eucD(CubeCoord a, CubeCoord b) {
	auto ap = euclideanCenter(a);
	auto bp = euclideanCenter(b);
	float d = (ap.x - bp.x) * (ap.x - bp.x) + (ap.y - bp.y) * (ap.y - bp.y);
	for (auto c : mirrorCenters) {
			CubeCoord a0 = cubeAdd(a, c);
			auto a0p = euclideanCenter(a0);
			d = min(d, (a0p.x - bp.x) * (a0p.x - bp.x) + (a0p.y - bp.y) * (a0p.y - bp.y));
		}
	return d;
}
void drawHex(CubeCoord c) {
	ofBeginShape();
	auto center = euclideanCenter(c);
	for (int theta = 30; theta < 360; theta += 60) {
		ofVertex(center.x + hexSizePix * cos(theta * pi / 180.)+pixelCenter.x, center.y + hexSizePix * sin(theta * pi / 180.) + pixelCenter.y);
	}
	ofEndShape();
	
}
void printCube(CubeCoord a) {
	cout << a.x << "," << a.y << "," << a.z;
}
CubeCoord recenter(CubeCoord a) {
	if (cubeDistance(a, center) > hexSize) {
		for (auto c : mirrorCenters) {
			if (cubeDistance(a, c) <= hexSize) {
				return cubeSub(a, c);
			}
		}
	}
	return a;
}
vector <CubeCoord> neighbors(CubeCoord a) {
	vector<CubeCoord> n;
	for (CubeCoord dir : cubeDirections) {
		n.push_back(recenter(cubeAdd(a, dir)));
	}
	return n;
}
void centroid() {
	cout << "centroid" << endl;
	int pieceCount = 0;
	for (auto o : pieceOrigins) {
		int newX = 0;
		int newY = 0;
		int newZ = 0;
		int newCount = 0;
		for (int i = -hexSize; i <= hexSize; ++i) {
			for (int j = -hexSize; j <= hexSize; ++j) {
				if (abs(i + j) > hexSize) {
					continue;
				}
				CubeCoord c = { i,j,-i - j };
				CubeCoord c0 = c;
				if (hexes[c].first != pieceCount) {
					continue;
				}
				int d = cubeDistance(o, c0);
				
				for (auto mirror : mirrorCenters) {
					int d1 = cubeDistance(o, cubeAdd(c0, mirror));
					if (d1 < d) {
						c = cubeAdd(c0, mirror);
						d = d1;
					}
				}
				
				newX += c.x;
				newY += c.y;
				newZ += c.z;
				newCount += 1;
			
			}
		}
		if (newCount == 0) {
			continue;
		}
		pieceOrigins[pieceCount] = recenter({ newX / newCount,newY / newCount,-newX / newCount - newY / newCount });
		pieceCount += 1;
	}
	cout << "done centroid" << endl;
}
void voronoi() {
	cout << "voronoi" << endl;
	for (int i = -hexSize; i <= hexSize; ++i) {
		for (int j = -hexSize; j <= hexSize; ++j) {
			if (abs(i + j) > hexSize) {
				continue;
			}
			CubeCoord c = { i,j,-i - j };
			hexes[c] = make_pair(-1, 100000000000);
		}
	}
	int pieceCount = 0;
	for (CubeCoord o : pieceOrigins) {
		hexes[o] = make_pair(pieceCount, 0);
		auto currentNeighbors = neighbors(o);
		while (currentNeighbors.size() > 0) {
			auto neighbor = currentNeighbors[currentNeighbors.size() - 1];
			currentNeighbors.pop_back();
			//int d = squareCubeEuclideanDistance(neighbor, o);
			float d = eucD(neighbor, o);
			auto check = hexes[neighbor];
			if (check.second > d && check.first != pieceCount) {
				hexes[neighbor] = make_pair(pieceCount, d);
				auto newNeighbors = neighbors(neighbor);
				for (auto newNeighbor : newNeighbors) {
					if (hexes[newNeighbor].first != pieceCount) {
						currentNeighbors.push_back(newNeighbor);
					}
				}
			}
		}
		pieceCount += 1;
	}
	cout << "done voronoi" << endl;
}
void testGrid() {
	cout << "testgrid" << endl;
	for (int i = 0; i < 200; ++i) {
		int x = ofRandom(-hexSize, hexSize);
		int y = ofRandom(-hexSize, hexSize);
		CubeCoord c = { x,y,-x - y };
		if (abs(c.z) > hexSize) {
			continue;
		}
		bool good = true;
		for (auto o : pieceOrigins) {
			if (squareCubeEuclideanDistance(o, c) < 10) {
				good = false;
				break;
			}
		}
		if (good) {
			pieceOrigins.push_back(c);
		}
	}
	
	//for(int i=0;i<1000;++i){
	//voronoi();
	//centroid();
	//}
	
	voronoi();
	cout << "done testgrid" << endl;
}

ofFbo test;
void drawImage() {
	cout << "start make image";
	test.allocate(4000, 4000);
	test.begin();
	for (int i = -hexSize; i <= hexSize; ++i) {
		for (int j = -hexSize; j <= hexSize; ++j) {
			if (abs(i + j) > hexSize) {
				continue;
			}
			CubeCoord c = { i,j,-i - j };
			auto h = hexes[c].first;
			if (h == -1) {
				cout << "didn't find one";
				ofSetColor({ 0,0,0,0 });
			}
			else {
				ofSetColor(colors[h]);
			}
			drawHex(c);
		}
	}
	
	for(auto mirror:mirrorCenters){
	for (int i = -hexSize; i <= hexSize; ++i) {
		for (int j = -hexSize; j <= hexSize; ++j) {
			if (abs(i + j) > hexSize) {
				continue;
			}
			CubeCoord c = { i,j,-i - j };
			auto h = hexes[c].first;
			if (h == -1) {
				cout << "didn't find one";
				ofSetColor({ 0,0,0,0 });
			}
			else{
			ofSetColor(colors[h]);
			}
			drawHex(cubeAdd(c,mirror));
		}
	}
	}
	
	test.end();
	cout << "done make image";
}

//--------------------------------------------------------------
void ofApp::setup(){
	CubeCoord mirror1 = { 2 * hexSize + 1,-hexSize - 1,-hexSize };
	mirrorCenters.push_back(mirror1);
	cout << "mirror 0 " << mirror1.x << "," << mirror1.y << "," << mirror1.z << endl;
	cout << "distance to mirror " << cubeDistance(center, mirror1) << endl;
	for (int i = 0; i < 5; ++i) {
		mirror1 = rotateClockwise(mirror1);
		cout << "mirror " << i + 1 << " " << mirror1.x << "," << mirror1.y << "," << mirror1.z << endl;
		mirrorCenters.push_back(mirror1);
	}
	CubeCoord b = { 2,0,-2 };
	CubeCoord c = { 0,-2,2 };
	cout << "cube distance " << cubeDistance(b, c) << endl;
	cout << "sq distance tiled" << squareCubeEuclideanDistance(b, c) << endl;
	CubeCoord d = { 3,0,-3 };
	CubeCoord mirr = { 3,2,-5 };
	cout << "distance to mirror " << cubeDistance(d, mirr) << endl;
	cout << "recenter ";
	printCube(recenter(d));
	cout << endl;
	for (int i = 0; i < 1000; ++i) {
		ofColor c = { ofRandom(0,255),ofRandom(0,255) ,ofRandom(0,255) };
		colors.push_back(c);
	}
	testGrid();
	drawImage();
	
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){
	voronoi();
	centroid();
	drawImage();
	test.draw(0, 0);
	
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
