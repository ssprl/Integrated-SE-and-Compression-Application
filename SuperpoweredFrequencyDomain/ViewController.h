#import <UIKit/UIKit.h>

@interface ViewController: UIViewController{
@public
    float x;
    //,x1,x2,x3;
}
@property (weak,nonatomic) IBOutlet UILabel *betaLabel;
@property (weak,nonatomic) IBOutlet UISlider *betaSlider;
@property (weak, nonatomic) IBOutlet UILabel *snrLabel;
- (IBAction)betaValue:(id)sender;
- (IBAction)snrDisplay:(id)sender;
@property (strong, nonatomic) IBOutlet UIImageView *ssprlImage;
@property (strong, nonatomic) IBOutlet UISwitch *EnhancedSwitch;
- (IBAction)switchPressed:(id)sender;
@property (strong, nonatomic) IBOutlet UILabel *OnDisplay;
@property (weak, nonatomic) IBOutlet UILabel *musicalLabel;
@property (weak, nonatomic) IBOutlet UILabel *muuLabel;
@property (weak, nonatomic) IBOutlet UILabel *vLabel;
@property (weak, nonatomic) IBOutlet UISlider *muuSlider;
@property (weak, nonatomic) IBOutlet UISlider *vSlider;
@property (weak, nonatomic) IBOutlet UIStepper *musicalStepper;
- (IBAction)musicalValue:(id)sender;
@property (weak, nonatomic) IBOutlet UILabel *muu;
@property (weak, nonatomic) IBOutlet UILabel *v;

@property (strong, nonatomic) IBOutlet UIStepper *stepper;
- (IBAction)stepbeta:(id)sender;


@property (weak, nonatomic) IBOutlet UILabel *musical;
- (IBAction)muuValue:(id)sender;
- (IBAction)vValue:(id)sender;


@end

