#pragma once
namespace fastMarch {
	class Cell {
		public:
		int x, y;
		double t;
		int ID;
		int status;

		Cell();
		Cell(int _x, int _y, float _t, int _ID);
		int getX() const;
		int getY() const;
		float getT() const;
		void setT(float _t);
		void setID(int _ID);
		int getID() const;
		int getStatus() const;
		void setStatus(int _status);

	};
	class CellComparator {
	public:
		int operator() (const std::reference_wrapper<Cell> c1, const std::reference_wrapper<Cell> c2);
	};
	class March {
	public:
		int curID;
		int numPts;
		std::vector<std::vector<Cell>> cells;
		March(int w, int h, std::vector<int> & mask, std::vector<std::pair<int, int>> & origins);
	};
	std::vector<std::pair<int, int>> centroid(March & m);
}