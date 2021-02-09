#include <vector>
#include <queue>
#include "fastMarch.h"
#include <map>
#include <iostream>
using namespace std;
namespace fastMarch {

		Cell::Cell() {
			x = 0; y = 0;
			t = 1e9;
			ID = -2;
			status = 0;
		}
		Cell::Cell(int _x, int _y, float _t, int _ID) {
			x = _x;
			y = _y;
			t = _t;
			ID = _ID;
			status = 0;
		}
		int Cell::getX() const { return x; }
		int Cell::getY() const { return y; }
		float Cell::getT() const { return t; }
		void Cell::setT(float _t) { t = _t; }
		void Cell::setID(int _ID) { ID = _ID; }
		int Cell::getID() const { return ID; }
		int Cell::getStatus() const { return status; }
		void Cell::setStatus(int _status) { status = _status; }

		int CellComparator::operator() (const reference_wrapper<Cell> c1, const reference_wrapper<Cell> c2) {

			return remove_reference<Cell&>::type(c1).getT() > remove_reference<Cell&>::type(c2).getT();
		}



		

		March::March(int w, int h, vector<int> & mask, vector<pair<int, int>> & origins)
		{
			numPts = origins.size();
			curID = 0;
			const pair<int, int> neighbor[4] = { make_pair(-1,0),make_pair(1,0),make_pair(0,-1),make_pair(0,1) };
			cells.resize(w);
			for (int i = 0; i < w; i++) {
				cells[i].resize(h);
			}
			for (int x = 0; x < w; x++) {
				for (int y = 0; y < h; y++) {
					cells[x][y] = Cell(x, y, 1e9, -1);
				}
			}
			for (int i = 0; i < w; i++) {
				for (int j = 0; j < h; j++) {
					if (mask[j*w + i] >= 0) {
						cells[i][j].setID(-2);
					}
					cells[i][j].setStatus(0);
				}
			}
			int count = 0;
			for (pair<int, int> o : origins) {
				count += 2;
				int x = o.first;
				int y = o.second;
				//cout << "processing origin " << x << "," << y << endl;
				//vector <Cell> accepted;
				priority_queue <reference_wrapper<Cell>, vector<reference_wrapper<Cell>>, CellComparator> considered;
				//accepted.push_back(cells[x][y]);

				cells[x][y].setID(curID);
				cells[x][y].setT(0);
				for (pair<int, int> n : neighbor) {
					int xd = n.first + x;
					int yd = n.second + y;
					if (xd < w && xd >= 0 && yd < h && yd >= 0) {

						if (cells[xd][yd].getID() != -2) {
							cells[xd][yd].setT(1);
							cells[xd][yd].setID(curID);
							cells[xd][yd].setStatus(count);
							considered.push(cells[xd][yd]);
						}
					}
				}

				float Tx, Ty, a, b, c;
				float arg, Ttemp;
				while (!considered.empty()) {
					Cell acceptedCell = considered.top();
					float T = acceptedCell.getT();
					acceptedCell.setStatus(count+1);
					considered.pop();

					int x = acceptedCell.getX(); int y = acceptedCell.getY();
					//for neighbors of accepted cell
					for (pair<int, int> n : neighbor) {
						int cX = x + n.first;
						int cY = y + n.second;
						Tx = 0., Ty = 0.;
						//cout << cX << "," << cY << endl;

						if (cX < w && cX >= 0 && cY < h && cY >= 0) {
							Cell& look = cells[cX][cY];

							if (look.getT() > T && look.getID() != -2 && look.getStatus() < count) {
								//cout << "look T " << look.getT() << " accepted T " << T << endl;
								//cout << "setting id " << curID << endl;
								look.setID(curID);
								if (cX > 0) {

									Cell& lLeft = cells[cX - 1][cY];

									if (lLeft.getID() == curID) {
										Tx = lLeft.getT();

									}
								}
								if (cX < w - 1) {
									Cell& lRight = cells[cX + 1][cY];
									if (lRight.getID() == curID) {
										if (Tx > 0) {
											Tx = min(Tx, lRight.getT());
										}
										else {
											Tx = lRight.getT();
										}

									}
								}
								if (cY > 0) {
									Cell& lDown = cells[cX][cY - 1];
									if (lDown.getID() == curID) {
										Ty = lDown.getT();
									}
								}
								if (cY < h - 1) {
									Cell& lUp = cells[cX][cY + 1];
									if (lUp.getID() == curID) {
										if (Ty > 0) {
											Ty = min(Ty, lUp.getT());
										}
										else {
											Ty = lUp.getT();
										}
									}
								}
								//cout << "Delta X,Y " << Tx << " " << Ty << endl;

								b = (Tx + Ty) / 2;
								c = 2 * Tx * Ty - Tx * Tx - Ty * Ty + 2;

								Ttemp = 0.;

								if (c >= 0) {
									Ttemp = b + .5 * sqrt(c);
									//cout << Ttemp;
								}
								if (Ttemp > max(Tx, Ty)) {
									//cout << "setting T normal " << Ttemp << endl;
									look.setT(Ttemp);
								}
								else {
									if (Tx == 0) {
										Tx = 1e9;
									}
									if (Ty == 0) {
										Ty = 1e9;
									}
									//cout << "SETTING T fallback " << min(Tx, Ty) + 1 <<endl;
									look.setT(min(Tx, Ty) + 1);
								}
								look.setStatus(count);
								considered.push(look);
							}
						}

					}

				}
				curID += 1;
			}

		}

		std::vector<std::pair<int, int>> centroid(March & m)
		{
			int w = m.cells.size();
			int h = m.cells[0].size();
			vector <pair<int, int> > idLocations(m.numPts,make_pair(0,0));
			vector <int> idCount(m.numPts,0);
			
			for (int i = 0; i < w; i++) {
				for (int j = 0; j < h; j++) {
					int id = m.cells[i][j].getID();
					if (id < 0) {
						continue;
					}
					idCount[id] += 1;
					pair<int, int> & loc = idLocations[id];
					loc.first += i;
					loc.second += j;
					
					
				}
			}
			vector <pair <int, int>> origins;
			for (int i=0;i<m.numPts;++i) {
				int count = idCount[i];
				pair <int, int> & loc = idLocations[i];
				if (count > 0) {
					loc.first = int(loc.first / count);
					loc.second = int(loc.second / count);
					if (m.cells[loc.first][loc.second].getID() == -2) {
						cout << "bad centroid: find closest" << endl;
						int minD = 1e10;
						int minX = loc.first;
						int minY = loc.second;
						for (int x = 0; x < w; x++) {
							for (int y = 0; y < h; y++) {
								int id = m.cells[x][y].getID();
								if (id == i) {
									int d = (x - loc.first)*(x - loc.first) + (y - loc.second)*(y - loc.second);
									if (d < minD) {
										minD = d;
										minX = x;
										minY = y;
									}
								}
							}
						}
						loc.first = minX;
						loc.second = minY;
					}
				}
				origins.push_back(loc);

			}
			return origins;
		}

}
