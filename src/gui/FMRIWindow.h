/////////////////////////////////////////////////////////////////////////////
// Name:            FMRIWindow.h
// Author:          ---
// Creation Date:   ---
//
// Description: This is the implementation file for the RT-fMRI exploration method.
/////////////////////////////////////////////////////////////////////////////


#ifndef FMRIWINDOW_H_
#define FMRIWINDOW_H_

#include "MainCanvas.h"
#include "MyListCtrl.h"

#include "../misc/Algorithms/Helper.h"

#include <wx/scrolwin.h>
#include <wx/statline.h>

class MainFrame;
class wxToggleButton;

class FMRIWindow: public wxScrolledWindow
{
public:
    FMRIWindow(){};
    FMRIWindow( wxWindow *pParent, MainFrame *pMf, wxWindowID id, const wxPoint &pos, const wxSize &size );

    ~FMRIWindow(){};
    void OnPaint( wxPaintEvent &event );
    void OnSize( wxSizeEvent &event );
    wxSizer* getWindowSizer();
	void SetSelectButton(); 
	void SetStartButton(){ m_pBtnStart->Enable( true ); m_pBtnStart->SetBackgroundColour(wxColour( 147, 255, 239 ));}

public:

	void OnStartRTFMRI					       ( wxCommandEvent& event );
	void OnSliderCorrThreshMoved			   ( wxCommandEvent& event );
	void OnSliderSizePMoved					   ( wxCommandEvent& event );
	void OnSliderAlphaMoved					   ( wxCommandEvent& event );
	void onConvertRestingState				   ( wxCommandEvent& event );
    void OnSelectClusters                      ( wxCommandEvent& event );
	
	
	
private:
    MainFrame           *m_pMainFrame;
	wxButton            *m_pBtnSelectFMRI;
	wxToggleButton      *m_pBtnStart;
	wxSlider            *m_pSliderCorrThreshold;
	wxTextCtrl          *m_pTxtCorrThreshBox;
    wxStaticText        *m_pTextCorrThreshold;
	wxSlider            *m_pSliderSizeP;
	wxTextCtrl          *m_pTxtSizePBox;
    wxStaticText        *m_pTextSizeP;
	wxSlider            *m_pSliderAlpha;
	wxTextCtrl          *m_pTxtAlphaBox;
    wxStaticText        *m_pTextAlpha;
	wxButton			*m_pBtnConvertFMRI;
    wxButton            *m_pBtnSelectClusters;

	bool showRawData;
    
private:
    wxSizer *m_pFMRISizer;
    FMRIWindow( wxWindow *pParent, wxWindowID id, const wxPoint &pos, const wxSize &size );
    DECLARE_DYNAMIC_CLASS( FMRIWindow )
    DECLARE_EVENT_TABLE()
};

#endif /*FMRIWINDOW*/
