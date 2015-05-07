#ifndef TAGDETECTOR_H
#define TAGDETECTOR_H

#include <vector>
#include <memory>

#include "opencv2/opencv.hpp"

#include "apriltags/TagDetection.h"
#include "apriltags/TagFamily.h"

namespace apriltags {

	// Temporary
	struct TagDetectorParameters
	{
		// Quad params
		unsigned int minClusterPixels;
		unsigned int maxNumMaxima;
		float criticalAngle; // In radians
		float maxLineFitMSE;
		int minWhiteBlackDiff;
		bool deglitch;
		
		// Other params
		float quadDecimate;
		float quadSigma;
		bool refineEdges;
		bool refineDecode;
		bool refinePose;
		
		TagDetectorParameters() 
			: minClusterPixels( 5 ), maxNumMaxima( 10 ), criticalAngle( M_PI/18 ),
			maxLineFitMSE( 1.0 ), minWhiteBlackDiff( 15 ), deglitch( false ),
			refineEdges( true ), refinePose( false ), refineDecode( false )
		{}
	};
	
    class TagDetector {
    public:

        typedef std::shared_ptr<TagDetector> Ptr;
        
        const TagFamily thisTagFamily;

        //! Constructor
        // note: TagFamily is instantiated here from TagCodes
        TagDetector(const TagCodes& tagCodes) 
			: thisTagFamily(tagCodes) 
		{}

        std::vector<TagDetection> ExtractTags(const cv::Mat& image) const;

    };
	
	/*! \brief Class that efficiently performs multi-family tag detections. */
	class MultiTagDetector {
	public:
		
		typedef std::shared_ptr<MultiTagDetector> Ptr;

		const std::vector<TagFamily> tagFamilies;
		
		MultiTagDetector( const std::vector<TagFamily>& families ) : tagFamilies( families ) {}
		
		std::vector<TagDetection> extractTags( const cv::Mat& image ) const;
		
	};

} // namespace

#endif
