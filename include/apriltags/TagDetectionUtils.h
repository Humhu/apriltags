#ifndef TAG_DETECTION_UTILS_H

#include "apriltags/TagDetection.h"
#include "apriltags/TagFamily.h"
#include "apriltags/FloatImage.h"
#include "apriltags/Quad.h"
#include "apriltags/UnionFindSimple.h"
#include "apriltags/Segment.h"

namespace AprilTags
{
	
		
	/*! \brief Steps 1-2. Calculate image gradients */
	void preprocess_image( const FloatImage& orig, FloatImage& fim,
						   FloatImage& fimSeg, FloatImage& fimTheta, FloatImage& fimMag,
						   float sigma = 0.0, float segSigma = 0.8 );
		
	/*! \brief Step 3. Find edges */
	void extract_edges( const FloatImage& fimSeg, const FloatImage& fimMag,
						const FloatImage& fimTheta, UnionFindSimple& uf );
	
	
	/*! \brief Step 4-5. Fit line segments */
	void fit_segments( const FloatImage& fimSeg, const FloatImage& fimMag,
					   const FloatImage& fimTheta, UnionFindSimple& uf,
					   std::vector<Segment>& segments );
	
	/*! \brief Step 6-7. Find quads */
	void find_quads( std::vector<Segment>& segments,
					 const FloatImage& fimOrig,
					 const std::pair<int,int>& opticalCenter,
					 std::vector<Quad>& quads );
	
	/*! \brief Step 8-9. Decode quads from a float image for a given family. */
	std::vector<TagDetection> decode_quads( const FloatImage& fim,
											const std::vector<Quad>& quads,
											const TagFamily& family );
}

#endif
