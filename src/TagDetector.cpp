#include "apriltags/TagDetector.h"

namespace apriltags
{
	
std::vector<TagDetection> TagDetector::ExtractTags( const cv::Mat& im_orig ) const
{
	
	///////////////////////////////////////////////////////////
    // Step 1. Detect quads according to requested image decimation
    // and blurring parameters.
	cv::Mat quad_im;
	if( quadDecimate > 1.0 )
	{
		cv::resize( im_orig, quad_im, cv::Size(), 1.0/quadDecimate, 1.0/quadDecimate,
			cv::INTER_LINEAR );
	}
	
	if( quadSigma != 0 )
	{
		// compute a reasonable kernel width by figuring that the
        // kernel should go out 2 std devs.
        //
        // max sigma          ksz
        // 0.499              1  (disabled)
        // 0.999              3
        // 1.499              5
        // 1.999              7
		
		// TODO Look at these parameters and simplify
		float sigma = std::abs( quadSigma );
		int ksz = std::round( 4*sigma );
		cv::Size ksize( ksz, ksz );

		if( (ksz & 1) == 0 )
		{
			if( quadSigma > 0 )
			{
				// Apply a blur 
				cv::GaussianBlur( quad_im, quad_im, ksize, sigma );
				gaussian_blur( quad_im, sigma, ksz );
			}
			else
			{
				// Sharpen the image by sutracting the low frequency components
				cv::Mat blurred;
				cv::GaussianBlur( quad_im, blurred, ksize, sigma );
				quadSigma = quadSigma - blurred;
			}
		}
	}
	
	// Get quads
	// TODO
	std::vector<Quad> quads = QuadThresh( quad_im );
	
	// adjust centers of pixels so that they correspond to the
    // original full-resolution image.
	// TODO Move this into QuadThresh?
	if( quadDecimate > 1.0 )
	{
		BOOST_FOREACH( Quad& quad, quads )
		{
			for( unsigned int i = 0; i < 4; i++ )
			{
				quad.p[i].x *= quadDecimate;
				quad.p[i].y *= quadDecimate;
			}
		}
	}
	
	////////////////////////////////////////////////////////////////
    // Step 2. Decode tags from each quad.
    std::vector<TagDetection> detections;
	BOOST_FOREACH( const Quad& quad, quads )
	{
		// TODO
		TagDetection detection;
		if( DecodeQuad( im_orig, quad, detection ) )
		{
			detections.push_back( detection );
		}
	}
	
	////////////////////////////////////////////////////////////////
    // Step 3. Reconcile detections--- don't report the same tag more
    // than once. (Allow non-overlapping duplicate detections.)
    Polygon poly0, poly1;
	
	TagDetections goodDetections;
	char status[ detections.size() ];
	memset( status, 1, detections.size()*sizeof(char) );
	
    for( unsigned int i = 0; i < detections.size(); i++ )
	{
		if( !status[i] ) { continue; }
		TagDetection& det0 = detections[0];
		
		poly0 = Polygon( det0.p );
		
		bool clear = true;
		for( unsigned int j = i+1; j < detections.size(); j++ )
		{
			if( !status[j] ) { continue; }
			TagDetection& det1 = detections[1];
			
			if( det0.id != det1.id || det0.family != det1.family ) { continue; }
			
			poly1 = Polygon( det1.p );
			if( poly0.Overaps( poly1 ) )
			{
				if( det0.hamming < det1.hamming ||
					det0.hamming == det1.hamming && det0.goodness > det1.goodness )
				{
					status[j] = 0;
				}
				else
				{
					status[i] = 0; // Not actually useful
					clear = false;
					break;
				}
			}
		}
		
		if( clear )
		{
			goodDetections.push_back( det0 );
		}
	}
	
    return goodDetections;
    
}
