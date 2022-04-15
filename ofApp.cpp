#include <string>
#include "ofApp.h"
#include "ofxSvg.h"
struct hexInfo {
	int squaredDistance;
	int pieceID;
};
struct CubeCoord {
	int x{0};
	int y{0};
	int z{0};

};
bool operator==(const CubeCoord& lhs, const CubeCoord& rhs)
{
	return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}
bool operator<(const CubeCoord& l, const CubeCoord& r) {
	return (l.x < r.x || (l.x == r.x && l.y < r.y));
}
struct Point {
	float x;
	float y;
};
const float pi = 3.14159265359;
//radius of the major hexagon
int hexSize = 24;

//minor hex size in pixels (center to vertex)
float hexSizePix = 6.;
float hexWidth = sqrt(3) * hexSizePix;
float hexHeight = 2 * hexSizePix;

//cube to cartesian coordinate vectors
float xCube[2] = { hexWidth,0};
float yCube[2] = { hexWidth / 2.,3. / 4. * hexHeight };

//pixel coordinate of center of major hexagon
ofPoint pixelCenter((hexSize+.5)* hexWidth, (hexSize+.5)* hexHeight);
CubeCoord center = { 0,0,0 };
//cardinal hex directions
const CubeCoord cubeDirections[6] = { {1,-1,0},{1,0,-1},{0,1,-1},{-1,1,0},{-1,0,1},{0,-1,1} };

//holds the centers of the adjacent major hexagons
vector<CubeCoord> mirrorCenters;
//hexagons to pieces and squared distances
map<CubeCoord, pair<int, int>> hexes;
//piece origins
vector<CubeCoord> pieceOrigins;
//number of hexes in each blob
vector<int> blobSize;

vector<ofColor> colors;
float temp = 100.;
float cooling = .99;
//total hexes in grid
int totalHexes = 0;
//total number of blobs
int totalBlobs = 0;

//Puzzle Constants
float puzzleWidthInches = 5;
float puzzleHeightInches = 5;
float pixelToPPI = 1.0; //PDF PPI = 72.0
float inchToPixel = 72.; 
vector<ofPolyline> whimsies;
vector<vector<CubeCoord>> whimsyHexes;
vector<CubeCoord> mousehex;
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
int mirrorDistance(CubeCoord a, CubeCoord b) {
	int d = cubeDistance(a,b);
	if (d == 1) return d;
	for (int i = 0; i < 6; ++i) {
		d = min(d, cubeDistance(a, cubeAdd(b, mirrorCenters[i])));
	}
	return d;
}
CubeCoord closestMirror(CubeCoord ref, CubeCoord other){
	int d = cubeDistance(ref, other);
	CubeCoord closest = other;
	for (int i = 0; i < 6; ++i) {
		CubeCoord added = cubeAdd(other, mirrorCenters[i]);
		int d1 = cubeDistance(ref, added);
		if (d1 < d) {
			closest = added;
			d = d1;
		}
	}
	return closest;
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
		ofVertex(center.x + hexSizePix * cos(theta * pi / 180.)+pixelCenter.x + 200, center.y + hexSizePix * sin(theta * pi / 180.) + pixelCenter.y + 200);
	}
	ofEndShape();
	
}
CubeCoord cubeRound(float q, float r, float s) {
	float qr = round(q);
	float rr = round(r);
	float sr = round(s);

	float q_diff = abs(q - qr);
	float r_diff = abs(r - rr);
	float s_diff = abs(s - sr);

	if (q_diff > r_diff and q_diff > s_diff) {
		qr = -rr - sr;
	}
	else if (r_diff > s_diff) {
		rr = -qr - sr;
	}
	else
		sr = -qr - rr;
	CubeCoord c;
	c.x = int(qr); c.y = int(rr); c.z = int(sr);
	return c;
}
CubeCoord pixelToHex(ofPoint p){
	p = p - pixelCenter;
	float q = (sqrt(3.) / 3. * p.x - 1. / 3 * p.y) / hexSizePix;
	float r = (-2. / 3 * p.y) / hexSizePix;
	return cubeRound(q, r, -q - r);
}
vector<CubeCoord> whimsyIntersects(ofPolyline whim) {
	vector<CubeCoord> v;
	ofFbo offscreen;
	offscreen.allocate(puzzleWidthInches * inchToPixel, puzzleHeightInches * inchToPixel);
	ofSetColor(0);
	offscreen.begin();
	ofBeginShape();
	for (auto vert : whim) {
		ofVertex(vert);
	}
	ofEndShape();
	ofImage im;
	offscreen.readToPixels(im);
	offscreen.end();
	im.update();
	im.save("test.png");
	ofColor c;
	for (int i = 0; i < im.getWidth(); ++i) {
		for (int j = 0; j < im.getHeight(); ++j) {
			c = im.getColor(i, j);
			if (c == ofColor(0)) {
				CubeCoord coord = pixelToHex(ofPoint(i, j));
				if (find(v.begin(), v.end(), coord) == v.end()) {
					v.push_back(coord);
					cout << "hit " << coord.x << "," << coord.y << "," << coord.z << endl;
				}
			}
		}
	}
	return v;
}
void drawHexEdge(pair<CubeCoord, CubeCoord> edge) {
	CubeCoord a = edge.first;
	CubeCoord b = edge.second;
	if (cubeDistance(a, b) != 1) {
		cout << "bad edge";
		return;
	}
	auto center = euclideanCenter(a);
	//stupid way of doing this
	vector<ofPoint> aPoint;
	vector<ofPoint> bPoint;
	for (int theta = 30; theta < 360; theta += 60) {
		aPoint.push_back(ofPoint(center.x + hexSizePix * cos(theta * pi / 180.) + pixelCenter.x, center.y + hexSizePix * sin(theta * pi / 180.) + pixelCenter.y));
	}
	center = euclideanCenter(b);
	for (int theta = 30; theta < 360; theta += 60) {
		bPoint.push_back(ofPoint(center.x + hexSizePix * cos(theta * pi / 180.) + pixelCenter.x, center.y + hexSizePix * sin(theta * pi / 180.) + pixelCenter.y));
	}
	vector<ofPoint> edgePoints; //woooow
	for (int i = 0; i < 6; ++i) {
		ofPoint a = aPoint[i];
		for (int j = 0; j < 6; ++j) {
			ofPoint b = bPoint[j];
			if (a.distance(b) < .01) {
				edgePoints.push_back(a);
			}
		}
	}
	if (edgePoints.size() != 2) cout << "what are you doing?" << endl;
	ofDrawLine(edgePoints[0]+200, edgePoints[1]+200);
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
vector <CubeCoord> neighbors(CubeCoord a,bool noWhimsyNeighbor = false) {
	vector<CubeCoord> n;
	for (CubeCoord dir : cubeDirections) {
		CubeCoord c1 = recenter(cubeAdd(a, dir));
		if (noWhimsyNeighbor) {
			bool good = true;
			for (auto whim : whimsyHexes) {
				if (find(whim.begin(), whim.end(), c1) != whim.end()) {
					good = false;
					break;
				}
			}
			if (good) {
				n.push_back(c1);
			}

		}
		else {
			n.push_back(c1);
		}
	}
	return n;
}
vector<CubeCoord> unrecenterNeighbors(CubeCoord a) {
	vector<CubeCoord> n;
	for (CubeCoord dir : cubeDirections) {
		n.push_back(cubeAdd(a, dir));
	}
	return n;
}
void centroid() {
	//cout << "centroid" << endl;
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
	//cout << "done centroid" << endl;
}
void voronoi() {
	cout << "voronoi" << endl;
	for (int i = -hexSize; i <= hexSize; ++i) {
		for (int j = -hexSize; j <= hexSize; ++j) {
			if (abs(i + j) > hexSize) {
				continue;
			}
			CubeCoord c = { i,j,-i - j };
			bool inWhimsy = false;
			int whimsyID = -1;
			for(int i=0;i<whimsyHexes.size();i++){
				auto whim = whimsyHexes[i];
				if (find(whim.begin(), whim.end(), c) != whim.end()) {
					whimsyID = i;
					inWhimsy = true;
					break;
				}
			}
			int d = 100000000000;
			if (inWhimsy) d = -1;
			hexes[c] = make_pair(whimsyID, d);
		}
	}
	int pieceCount = 0;
	for (CubeCoord o : pieceOrigins) {
		hexes[o] = make_pair(pieceCount, 0);
		auto currentNeighbors = neighbors(o,true);
		while (currentNeighbors.size() > 0) {
			auto neighbor = currentNeighbors[currentNeighbors.size() - 1];
			currentNeighbors.pop_back();
			//int d = squareCubeEuclideanDistance(neighbor, o);
			float d = eucD(neighbor, o);
			auto check = hexes[neighbor];
			if (check.second > d && check.first != pieceCount) {
				hexes[neighbor] = make_pair(pieceCount, d);
				auto newNeighbors = neighbors(neighbor,true);
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

void countBlobs() {
	blobSize.clear();
	blobSize.resize(pieceOrigins.size());
	for (int i = -hexSize; i < hexSize+1; ++i) {
		for (int j = -hexSize; j < hexSize+1; ++j) {
			if (abs(i + j) > hexSize) {
				continue;
			}
			CubeCoord c = { i,j,-i - j };
			int C = hexes[c].first;
			if (C > -1) {
				totalHexes += 1;
				blobSize[C] += 1;
			}
		}
	}
}

void testGrid() {
	cout << "testgrid" << endl;
	for (int i = 0; i < 200000; ++i) {
		int x = ofRandom(-hexSize, hexSize);
		int y = ofRandom(-hexSize, hexSize);
		CubeCoord c = { x,y,-x - y };
		if (abs(c.z) > hexSize) {
			continue;
		}
		bool good = true;
		for (auto o : pieceOrigins) {
			if (mirrorDistance(o, c) < 6) {
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
	cout << "num pieces: " <<  pieceOrigins.size() << endl;
}

ofFbo test;
void drawImage(bool save = false) {

	test.allocate(hexWidth*hexSize*3, hexHeight * hexSize * 3);
	if (save) {
		ofBeginSaveScreenAsPDF("test image.pdf");
	}
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
	if (save) {
		ofEndSaveScreenAsPDF();
	}
}
pair<CubeCoord, CubeCoord> makeHexPair(CubeCoord a, CubeCoord b) {
	if (a.x < b.x) return make_pair(a, b);
	if (a.x > b.x) return make_pair(b, a);
	if (a.y < b.y) return make_pair(a, b);
	return make_pair(b, a);
}
void exportPdf() {
	
	vector<pair<CubeCoord, CubeCoord>> drawnEdges;
	
	for (int i = 0; i < blobSize.size(); ++i) {
		vector<CubeCoord> blobHexes;
		blobHexes.clear();
		for (int j = -hexSize; j <= hexSize; ++j) {
			for (int k = -hexSize; k <= hexSize; ++k) {
				if (abs(j + k) > hexSize) continue;
				CubeCoord c = { j,k,-j - k };
				if (hexes[c].first == i) blobHexes.push_back(c);
			}
		}
		if (blobHexes.size() == 0) continue;
		CubeCoord pieceOrigin = blobHexes[0];
		for (int j = 0; j < blobHexes.size(); ++j) {
			blobHexes[j] = closestMirror(pieceOrigin, blobHexes[j]);
		}
		for (int j = 0; j < blobHexes.size(); ++j) {
			CubeCoord curr = blobHexes[j];
			vector<CubeCoord> neighs = unrecenterNeighbors(curr);
			for (auto neigh : neighs) {
				if (find(blobHexes.begin(), blobHexes.end(), neigh) == blobHexes.end()) {
					pair<CubeCoord, CubeCoord> hexEdge = makeHexPair(curr, neigh);
					if (find(drawnEdges.begin(), drawnEdges.end(), hexEdge) == drawnEdges.end()) {
						drawnEdges.push_back(hexEdge);
					}
				}
			}
		}
	}
	ofFbo offscreen;
	offscreen.allocate(hexSize * hexHeight * 4., hexSize * hexHeight * 4.,OF_IMAGE_COLOR);
	offscreen.begin();
	ofBeginSaveScreenAsPDF("test.pdf");
	drawHex(center);
	ofSetColor(255, 0, 255);
	for (auto hexEdge : drawnEdges) {
		drawHexEdge(hexEdge);
	}
	offscreen.end();
	ofEndSaveScreenAsPDF();
}

int countSame(CubeCoord a,vector<CubeCoord>& neighs, vector<CubeCoord>&others) {
	int count = 0;
	neighs.clear();
	others.clear();
	auto n = neighbors(a);
	int me = hexes[a].first;
	if (me == -1) return -1;
	for (auto neigh : n) {
		int them = hexes[neigh].first;
		if (them == me) {
			count += 1; neighs.push_back(neigh);
		}
		else {
			others.push_back(neigh);
		}

	}
	return count;
}
bool checkContinuous(vector<CubeCoord> neighbors) {

	if (neighbors.size() == 0) return true;
	int c0 = neighbors.size();
	CubeCoord neigh = neighbors[neighbors.size()-1];
	neighbors.pop_back();
	queue<CubeCoord> goodNeighbors;
	goodNeighbors.push(neigh);
	int countGood = 1;
	while(goodNeighbors.size()>0){
		CubeCoord cur = recenter(goodNeighbors.front());
		goodNeighbors.pop();
		for (int i = neighbors.size()-1; i >= 0;--i) {
			CubeCoord other = neighbors[i];
			if (mirrorDistance(cur, other) == 1) {
				countGood += 1;
				goodNeighbors.push(other);
				neighbors.erase(neighbors.begin() + i);
			}
		}
		
	}
	if (countGood == c0) return true;
	return false;
}
bool flip() {
	int i = ofRandom(-hexSize-1, hexSize+1);
	int j = ofRandom(-hexSize-1, hexSize+1);
	
	while (abs(i + j) > hexSize) {
		i = ofRandom(-hexSize - 1, hexSize + 1);
		j = ofRandom(-hexSize - 1, hexSize + 1);
	}
	int k = -i - j;
	CubeCoord c = { i,j,k };
	c = recenter(c);
	int meme = hexes[c].first;
	vector<CubeCoord> neighbors;
	vector<CubeCoord> others;
	
	int count = countSame(c, neighbors, others);
	if (neighbors.size() < 3) return false;
	if (neighbors.size() <= others.size()) {
		int d = others.size() - neighbors.size();
		if (d * 10 > temp+ofRandom(1)-.5) { cout << "temp" << endl; return false;
	}
	}

	
	if (neighbors.size() == 6) {
		//cout << "same neighbors" << endl;
		return false;
	}
	if (checkContinuous(neighbors)) {
		CubeCoord other = others[int(ofRandom(others.size()))];
		int newColor = hexes[other].first;
		if (blobSize[meme] < 6) return false;
		if (blobSize[meme]*1.5 < blobSize[newColor]) {
			return false;
		}
		neighbors.clear(); others.clear();
		if (countSame(other, neighbors, others) < 3) return false;
		if (!checkContinuous(neighbors)) return false;
		//if (mirrorDistance(pieceOrigins[newColor], c) > 4*2) return false;
		pair<int, int> thisOne = hexes[c];
		blobSize[meme] -= 1;
		blobSize[newColor] += 1;
		hexes[c] = make_pair(newColor, thisOne.second);
		//cout << "flipped. was " << meme << " now I'm " << newColor << endl;
		return true;
	}
	//cout << "noncontinous neighbors " << endl;
	return false;

}
//--------------------------------------------------------------
void ofApp::setup(){
	ofxSVG svg;
	svg.load("cactus.svg");
	ofPath path = svg.getPaths()[0];
	path.setStrokeWidth(1);
	ofPolyline pline = path.getOutline()[0];
	
	pline.scale(1, 1);

	pline.translate(ofPoint(inchToPixel * 1, inchToPixel * 1));
	
	whimsies.push_back(pline);
	whimsyHexes.push_back(whimsyIntersects(whimsies[0]));

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
	for (int i = 0; i < 10000; ++i) {
		ofColor c = { ofRandom(0,255),ofRandom(0,255) ,ofRandom(0,255) };
		colors.push_back(c);
	}
	testGrid();
	drawImage();
	
}
int countFrame = 0;
float iter = 100;
float increase = 1.01;
//--------------------------------------------------------------
void ofApp::update(){
	
	int iters = 10000;
	int countChange = 0;
	temp = temp * cooling;
	countFrame += 1;
	bool changed = 0;
	if (countFrame < 5) {
		voronoi();
		centroid();
		countBlobs();
	}

	else {
		centroid();
		for (int i = 0; i < iters; ++i) {
			countChange += int(flip());
			
		}
	}
	//printf("\33[2K\r");
	cout << temp << ":" << countChange*1.0/ iter << endl;
	if (countChange>0) {
		drawImage(countFrame%100==0);
	}
	test.begin();
	for (auto whim : whimsies) {
		ofSetColor(0, 0, 0);
		ofBeginShape();
		for (auto vert : whim) {
			ofVertex(vert+pixelCenter);
		}
		ofEndShape();
	}
	for (auto whim : whimsyHexes) {

		ofSetColor(0, 255, 0);

		for (auto hex : whim) {
			cout << hex.x << "," << hex.y << "," << hex.z << endl;
			drawHex(hex);
		}
	}
	ofSetColor(255, 0, 255);
	for (auto mouse : mousehex) {
		drawHex(mouse);
	}
	test.end();
	if ( countFrame%100==0) {
		ofImage i; i.allocate(test.getWidth(), test.getHeight(), OF_IMAGE_COLOR);
		test.readToPixels(i);
		i.update();
		//i.resize(i.getWidth() / 4., i.getHeight() / 4.);
		auto s = to_string(countFrame);
		auto pad = string(6 - s.size(), '0');
		i.save("test" +pad+s + ".png");
		exportPdf();
	}
}

//--------------------------------------------------------------
void ofApp::draw(){
	
	ofSetColor(255);
	test.draw(0, 0);
	
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
	if (key == 's') exportPdf();
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
	mousehex.push_back(pixelToHex(ofPoint(x, y)));
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
