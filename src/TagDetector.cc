#include <algorithm>
#include <cmath>
#include <climits>
#include <map>
#include <vector>
#include <iostream>

#include <Eigen/Dense>

#include "apriltags/TagDetector.h"
#include "apriltags/TagDetectionUtils.h"

#include <boost/date_time/posix_time/posix_time.hpp>

using namespace std;

namespace AprilTags {
	
std::vector<TagDetection> TagDetector::extractTags(const cv::Mat& image) const {

	typedef boost::posix_time::ptime Time;
    typedef boost::posix_time::time_duration TimeDuration;
	
	// convert to internal AprilTags image (todo: slow, change internally to OpenCV)
	int width = image.cols;
	int height = image.rows;
	
	std::pair<int,int> opticalCenter(width/2, height/2);
	cv::Mat fimOrigCV;
	
// 	Time start = boost::posix_time::microsec_clock::universal_time();
	image.convertTo( fimOrigCV, CV_32FC1, 1.0/255.0 );
// 	Time finish = boost::posix_time::microsec_clock::universal_time();
// 	TimeDuration convertTime = finish - start;
// 	std::cout << "Took " << convertTime << " to convert image." << std::endl;
	
	//================================================================
	// Step one: preprocess image (convert to grayscale) and low pass if necessary
	
	//================================================================
	// Step two: Compute the local gradient. We store the direction and magnitude.
	// This step is quite sensitve to noise, since a few bad theta estimates will
	// break up segments, causing us to miss Quads. It is useful to do a Gaussian
	// low pass on this step even if we don't want it for encoding.
	cv::Mat fimCV, fimSegCV, fimThetaCV, fimMagCV;
// 	start = boost::posix_time::microsec_clock::universal_time();
	preprocess_image( fimOrigCV, fimCV, fimSegCV, fimThetaCV, fimMagCV );
// 	finish = boost::posix_time::microsec_clock::universal_time();
// 	TimeDuration preprocessTime = finish - start;
// 	std::cout << "Took " << preprocessTime << " to process image." << std::endl;
	
	//================================================================
	// Step three. Extract edges by grouping pixels with similar
	// thetas together. This is a greedy algorithm: we start with
	// the most similar pixels.  We use 4-connectivity.
	UnionFindSimple uf( fimSegCV.cols*fimSegCV.rows );
// 	start = boost::posix_time::microsec_clock::universal_time();
	extract_edges( fimSegCV, fimMagCV, fimThetaCV,  uf );
// 	finish = boost::posix_time::microsec_clock::universal_time();
// 	TimeDuration extractTime = finish - start;
// 	std::cout << "Took " << extractTime << " to extract edges." << std::endl;
	
	// Steps 4-5
	std::vector<Segment> segments;
// 	start = boost::posix_time::microsec_clock::universal_time();
	fit_segments( fimSegCV, fimMagCV, fimThetaCV, uf, segments );
// 	finish = boost::posix_time::microsec_clock::universal_time();
// 	TimeDuration fitTime = finish - start;
// 	std::cout << "Took " << fitTime << " to fit segments." << std::endl;
	
	// Steps 6-7
	std::vector<Quad> quads;
// 	start = boost::posix_time::microsec_clock::universal_time();
	find_quads( segments, fimOrigCV.size(), opticalCenter, quads );
// 	finish = boost::posix_time::microsec_clock::universal_time();
// 	TimeDuration quadTime = finish - start;
// 	std::cout << "Took " << quadTime << " to find quads." << std::endl;
	
	//================================================================
	// Step eight. Decode the quads. For each quad, we first estimate a
	// threshold color to decide between 0 and 1. Then, we read off the
	// bits and see if they make sense.
	//================================================================
	//Step nine: Some quads may be detected more than once, due to
	//partial occlusion and our aggressive attempts to recover from
	//broken lines. When two quads (with the same id) overlap, we will
	//keep the one with the lowest error, and if the error is the same,
	//the one with the greatest observed perimeter.
// 	start = boost::posix_time::microsec_clock::universal_time();
	std::vector<TagDetection> dets = decode_quads( fimCV, quads, thisTagFamily );
// 	finish = boost::posix_time::microsec_clock::universal_time();
// 	TimeDuration decodeTime = finish - start;
// 	std::cout << "Took " << decodeTime << " to decode quads." << std::endl;
	
	return dets;
}

} // namespace
