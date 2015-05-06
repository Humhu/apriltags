#include "apriltags/TagDetectionUtils.h"

#include "apriltags/Gaussian.h"
#include "apriltags/Edge.h"
#include "apriltags/GLine2D.h"
#include "apriltags/GLineSegment2D.h"
#include "apriltags/GrayModel.h"
#include "apriltags/MathUtil.h"
#include "apriltags/Gridder.h"

namespace AprilTags
{
	
	void preprocess_image( const FloatImage& orig, FloatImage& fim, 
					   FloatImage& fimSeg, FloatImage& fimTheta, FloatImage& fimMag,
					   float sigma, float segSigma )
{
	fim = orig;

	//! Gaussian smoothing kernel applied to image (0 == no filter).
	/*! Used when sampling bits. Filtering is a good idea in cases
	* where A) a cheap camera is introducing artifical sharpening, B)
	* the bayer pattern is creating artifacts, C) the sensor is very
	* noisy and/or has hot/cold pixels. However, filtering makes it
	* harder to decode very small tags. Reasonable values are 0, or
	* [0.8, 1.5].
	*/

	//! Gaussian smoothing kernel applied to image (0 == no filter).
	/*! Used when detecting the outline of the box. It is almost always
	* useful to have some filtering, since the loss of small details
	* won't hurt. Recommended value = 0.8. The case where sigma ==
	* segsigma has been optimized to avoid a redundant filter
	* operation.
	*/

	if (sigma > 0) {
		int filtsz = ((int) max(3.0f, 3*sigma)) | 1;
		std::vector<float> filt = Gaussian::makeGaussianFilter(sigma, filtsz);
		fim.filterFactoredCentered(filt, filt);
	}	
	
	fimSeg = FloatImage();
	if (segSigma > 0) {
		if (segSigma == sigma) {
			fimSeg = fim;
		} else {
		// blur anew
		int filtsz = ((int) max(3.0f, 3*segSigma)) | 1;
		std::vector<float> filt = Gaussian::makeGaussianFilter(segSigma, filtsz);
		fimSeg = orig;
		fimSeg.filterFactoredCentered(filt, filt);
		}
	} else {
		fimSeg = orig;
	}

	fimTheta = FloatImage(fimSeg.getWidth(), fimSeg.getHeight());
	fimMag = FloatImage(fimSeg.getWidth(), fimSeg.getHeight());

	#pragma omp parallel for
	for (int y = 1; y < fimSeg.getHeight()-1; y++) {
		for (int x = 1; x < fimSeg.getWidth()-1; x++) {
		float Ix = fimSeg.get(x+1, y) - fimSeg.get(x-1, y);
		float Iy = fimSeg.get(x, y+1) - fimSeg.get(x, y-1);

		float mag = Ix*Ix + Iy*Iy;
	#if 0 // kaess: fast version, but maybe less accurate?
		// TODO Benchmark with libbams
		float theta = MathUtil::fast_atan2(Iy, Ix);
	#else
		float theta = atan2(Iy, Ix);
	#endif
		
		fimTheta.set(x, y, theta);
		fimMag.set(x, y, mag);
		}
	}
}
	
void extract_edges( const FloatImage& fimSeg, const FloatImage& fimMag,
					const FloatImage& fimTheta, UnionFindSimple& uf )
{

	int width = fimSeg.getWidth();
	int height = fimSeg.getHeight();
	
	std::vector<Edge> edges(width*height*4);
	size_t nEdges = 0;

// Bounds on the thetas assigned to this group. Note that because
// theta is periodic, these are defined such that the average
// value is contained *within* the interval.
	{ // limit scope of storage
		/* Previously all this was on the stack, but this is 1.2MB for 320x240 images
		* That's already a problem for OS X (default 512KB thread stack size),
		* could be a problem elsewhere for bigger images... so store on heap */
		std::vector<float> storage(width*height*4);  // do all the memory in one big block, exception safe
		float * tmin = &storage[width*height*0];
		float * tmax = &storage[width*height*1];
		float * mmin = &storage[width*height*2];
		float * mmax = &storage[width*height*3];
					
		for (int y = 0; y+1 < height; y++) {
		for (int x = 0; x+1 < width; x++) {
									
			float mag0 = fimMag.get(x,y);
			if (mag0 < Edge::minMag)
			continue;
			mmax[y*width+x] = mag0;
			mmin[y*width+x] = mag0;
									
			float theta0 = fimTheta.get(x,y);
			tmin[y*width+x] = theta0;
			tmax[y*width+x] = theta0;
									
			// Calculates then adds edges to 'vector<Edge> edges'
			Edge::calcEdges(theta0, x, y, fimTheta, fimMag, edges, nEdges);
									
			// TODO Would 8 connectivity help for rotated tags?
			// Probably not much, so long as input filtering hasn't been disabled.
		}
		}
					
		edges.resize(nEdges);
		std::stable_sort(edges.begin(), edges.end());
		Edge::mergeEdges(edges,uf,tmin,tmax,mmin,mmax);
	}
}

void fit_segments( const FloatImage& fimSeg, const FloatImage& fimMag, 
				   const FloatImage& fimTheta, UnionFindSimple& uf,
				   std::vector<Segment>& segments )
{
	//================================================================
	// Step four: Loop over the pixels again, collecting statistics for each cluster.
	// We will soon fit lines (segments) to these points.
	std::map<int, std::vector<XYWeight> > clusters;
	for (int y = 0; y+1 < fimSeg.getHeight(); y++) {
		for (int x = 0; x+1 < fimSeg.getWidth(); x++) {
		if (uf.getSetSize(y*fimSeg.getWidth()+x) < Segment::minimumSegmentSize)
		continue;

		int rep = (int) uf.getRepresentative(y*fimSeg.getWidth()+x);
		
		std::map<int, std::vector<XYWeight> >::iterator it = clusters.find(rep);
		if ( it == clusters.end() ) {
		clusters[rep] = std::vector<XYWeight>();
		it = clusters.find(rep);
		}
		std::vector<XYWeight> &points = it->second;
		points.push_back(XYWeight(x,y,fimMag.get(x,y)));
		}
	}

	//================================================================
	// Step five: Loop over the clusters, fitting lines (which we call Segments).
	std::map<int, std::vector<XYWeight> >::const_iterator clustersItr;
	for (clustersItr = clusters.begin(); clustersItr != clusters.end(); clustersItr++) {
		std::vector<XYWeight> points = clustersItr->second;
		GLineSegment2D gseg = GLineSegment2D::lsqFitXYW(points);

		// filter short lines
		float length = MathUtil::distance2D(gseg.getP0(), gseg.getP1());
		if (length < Segment::minimumLineLength)
		continue;

		Segment seg;
		float dy = gseg.getP1().second - gseg.getP0().second;
		float dx = gseg.getP1().first - gseg.getP0().first;

		float tmpTheta = std::atan2(dy,dx);

		seg.setTheta(tmpTheta);
		seg.setLength(length);

		// We add an extra semantic to segments: the vector
		// p1->p2 will have dark on the left, white on the right.
		// To do this, we'll look at every gradient and each one
		// will vote for which way they think the gradient should
		// go. This is way more retentive than necessary: we
		// could probably sample just one point!

		float flip = 0, noflip = 0;
		for (unsigned int i = 0; i < points.size(); i++) {
		XYWeight xyw = points[i];
		
		float theta = fimTheta.get((int) xyw.x, (int) xyw.y);
		float mag = fimMag.get((int) xyw.x, (int) xyw.y);

		// err *should* be +M_PI/2 for the correct winding, but if we
		// got the wrong winding, it'll be around -M_PI/2.
		float err = MathUtil::mod2pi(theta - seg.getTheta());

		if (err < 0)
		noflip += mag;
		else
		flip += mag;
		}

		if (flip > noflip) {
		float temp = seg.getTheta() + (float)M_PI;
		seg.setTheta(temp);
		}

		float dot = dx*std::cos(seg.getTheta()) + dy*std::sin(seg.getTheta());
		if (dot > 0) {
		seg.setX0(gseg.getP1().first); seg.setY0(gseg.getP1().second);
		seg.setX1(gseg.getP0().first); seg.setY1(gseg.getP0().second);
		}
		else {
		seg.setX0(gseg.getP0().first); seg.setY0(gseg.getP0().second);
		seg.setX1(gseg.getP1().first); seg.setY1(gseg.getP1().second);
		}

		segments.push_back(seg);
	}
}

void find_quads( std::vector<Segment>& segments, const FloatImage& fimOrig,
				 const std::pair<int,int>& opticalCenter, std::vector<Quad>& quads)
{
	
	int width = fimOrig.getWidth();
	int height = fimOrig.getHeight();
	
	// Step six: For each segment, find segments that begin where this segment ends.
	// (We will chain segments together next...) The gridder accelerates the search by
	// building (essentially) a 2D hash table.
	Gridder<Segment> gridder(0,0,width,height,10);

	// add every segment to the hash table according to the position of the segment's
	// first point. Remember that the first point has a specific meaning due to our
	// left-hand rule above.
	for (unsigned int i = 0; i < segments.size(); i++) {
		gridder.add(segments[i].getX0(), segments[i].getY0(), &segments[i]);
	}

	// Now, find child segments that begin where each parent segment ends.
	for (unsigned i = 0; i < segments.size(); i++) {
		Segment& parentseg = segments[i];
		
		//compute length of the line segment
		GLine2D parentLine(std::pair<float,float>(parentseg.getX0(), parentseg.getY0()),
				std::pair<float,float>(parentseg.getX1(), parentseg.getY1()));

		Gridder<Segment>::iterator iter = gridder.find(parentseg.getX1(), parentseg.getY1(), 0.5f*parentseg.getLength());
		while(iter.hasNext()) {
		Segment &child = iter.next();
		if (MathUtil::mod2pi(child.getTheta() - parentseg.getTheta()) > 0) {
		continue;
		}

		// compute intersection of points
		GLine2D childLine(std::pair<float,float>(child.getX0(), child.getY0()),
				std::pair<float,float>(child.getX1(), child.getY1()));

		std::pair<float,float> p = parentLine.intersectionWith(childLine);
		if (p.first == -1) {
		continue;
		}

		float parentDist = MathUtil::distance2D(p, std::pair<float,float>(parentseg.getX1(),parentseg.getY1()));
		float childDist = MathUtil::distance2D(p, std::pair<float,float>(child.getX0(),child.getY0()));

		if (max(parentDist,childDist) > parentseg.getLength()) {
		// cout << "intersection too far" << endl;
		continue;
		}

		// everything's OK, this child is a reasonable successor.
		parentseg.children.push_back(&child);
		}
	}

	//================================================================
	// Step seven: Search all connected segments to see if any form a loop of length 4.
	// Add those to the quads list.
	std::vector<Segment*> tmp(5);
	for (unsigned int i = 0; i < segments.size(); i++) {
		tmp[0] = &segments[i];
		Quad::search(fimOrig, tmp, segments[i], 0, quads, opticalCenter);
	}
}

std::vector<TagDetection> decode_quads( const FloatImage& fim,
										const std::vector<Quad>& quads,
										const TagFamily& family )
{
	//================================================================
	// Step eight. Decode the quads. For each quad, we first estimate a
	// threshold color to decide between 0 and 1. Then, we read off the
	// bits and see if they make sense.
	std::vector<TagDetection> detections;
	int width = fim.getWidth();
	int height = fim.getHeight();

	for (unsigned int qi = 0; qi < quads.size(); qi++ ) {
		const Quad &quad = quads[qi];

		// Find a threshold
		GrayModel blackModel, whiteModel;
		const int dd = 2 * family.blackBorder + family.dimension;

		for (int iy = -1; iy <= dd; iy++) 
		{
			float y = (iy + 0.5f) / dd;
			for (int ix = -1; ix <= dd; ix++) 
			{
				float x = (ix + 0.5f) / dd;
				std::pair<float,float> pxy = quad.interpolate01(x, y);
				int irx = (int) (pxy.first + 0.5);
				int iry = (int) (pxy.second + 0.5);
				if (irx < 0 || irx >= width || iry < 0 || iry >= height)
					continue;
				float v = fim.get(irx, iry);
				if (iy == -1 || iy == dd || ix == -1 || ix == dd)
					whiteModel.addObservation(x, y, v);
				else if (iy == 0 || iy == (dd-1) || ix == 0 || ix == (dd-1))
					blackModel.addObservation(x, y, v);
			}
		}

		bool bad = false;
		unsigned long long tagCode = 0;
		for ( int iy = family.dimension-1; iy >= 0; iy-- ) 
		{
			float y = (family.blackBorder + iy + 0.5f) / dd;
			for (int ix = 0; ix < family.dimension; ix++ ) 
			{
				float x = (family.blackBorder + ix + 0.5f) / dd;
				std::pair<float,float> pxy = quad.interpolate01(x, y);
				int irx = (int) (pxy.first + 0.5);
				int iry = (int) (pxy.second + 0.5);
				if (irx < 0 || irx >= width || iry < 0 || iry >= height) 
				{
					// cout << "*** bad:  irx=" << irx << "  iry=" << iry << endl;
					bad = true;
					continue;
				}
				float threshold = (blackModel.interpolate(x,y) + whiteModel.interpolate(x,y)) * 0.5f;
				float v = fim.get(irx, iry);
				tagCode = tagCode << 1;
				
				if ( v > threshold)
					tagCode |= 1;
				}
		}

		if( bad )
			continue;
		
		TagDetection thisTagDetection;
		family.decode(thisTagDetection, tagCode);

		// compute the homography (and rotate it appropriately)
		Homography33 hom = quad.homography;
		hom.compute();
		thisTagDetection.homography = hom.getH();
		thisTagDetection.hxy = hom.getCXY();

		float c = std::cos(thisTagDetection.rotation*(float)M_PI/2);
		float s = std::sin(thisTagDetection.rotation*(float)M_PI/2);
		Eigen::Matrix3d R;
		R.setZero();
		R(0,0) = R(1,1) = c;
		R(0,1) = -s;
		R(1,0) = s;
		R(2,2) = 1;
		Eigen::Matrix3d tmp;
		tmp = thisTagDetection.homography * R;
		thisTagDetection.homography = tmp;

		// Rotate points in detection according to decoded
		// orientation.  Thus the order of the points in the
		// detection object can be used to determine the
		// orientation of the target.
		std::pair<float,float> bottomLeft = thisTagDetection.interpolate(-1,-1);
		int bestRot = -1;
		float bestDist = FLT_MAX;
		for ( int i=0; i<4; i++ ) 
		{
			float const dist = AprilTags::MathUtil::distance2D(bottomLeft, quad.quadPoints[i]);
			if ( dist < bestDist ) {
				bestDist = dist;
				bestRot = i;
			}
		}

		for (int i=0; i< 4; i++)
			thisTagDetection.p[i] = quad.quadPoints[(i+bestRot) % 4];

		if (thisTagDetection.good) 
		{
			thisTagDetection.cxy = quad.interpolate01(0.5f, 0.5f);
			thisTagDetection.observedPerimeter = quad.observedPerimeter;
			detections.push_back(thisTagDetection);
		}
	}

	//================================================================
	//Step nine: Some quads may be detected more than once, due to
	//partial occlusion and our aggressive attempts to recover from
	//broken lines. When two quads (with the same id) overlap, we will
	//keep the one with the lowest error, and if the error is the same,
	//the one with the greatest observed perimeter.

	std::vector<TagDetection> goodDetections;

	// NOTE: allow multiple non-overlapping detections of the same target.

	for( std::vector<TagDetection>::const_iterator it = detections.begin();
		it != detections.end(); 
		it++ ) 
	{
		const TagDetection &thisTagDetection = *it;

		bool newFeature = true;

		for ( unsigned int odidx = 0; odidx < goodDetections.size(); odidx++) 
		{
			TagDetection &otherTagDetection = goodDetections[odidx];

			if ( thisTagDetection.id != otherTagDetection.id ||
			! thisTagDetection.overlapsTooMuch(otherTagDetection) )
				continue;

			// There's a conflict.  We must pick one to keep.
			newFeature = false;

			// This detection is worse than the previous one... just don't use it.
			if ( thisTagDetection.hammingDistance > otherTagDetection.hammingDistance )
				continue;

			// Otherwise, keep the new one if it either has strictly *lower* error, or greater perimeter.
			if ( thisTagDetection.hammingDistance < otherTagDetection.hammingDistance ||
			thisTagDetection.observedPerimeter > otherTagDetection.observedPerimeter )
			{
				goodDetections[odidx] = thisTagDetection;
			}
		}

		if ( newFeature )
			goodDetections.push_back(thisTagDetection);

	}

	return goodDetections;
}

}
